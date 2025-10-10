#include <ngstd.hpp>
#include <core/statushandler.hpp>

namespace ngstd
{
  static mutex m;
  void BaseStatusHandler::PushStatus (const std::string & str)
  {
    lock_guard lock(m);
    // Ng_PushStatus(str);
    ngcore::PushStatus(str);
  }

  void BaseStatusHandler::PopStatus ()
  {
    lock_guard lock(m);
    // Ng_PopStatus();
    ngcore::PopStatus();    
  }

  void BaseStatusHandler::SetThreadPercentage (double percent)
  {
    // Ng_SetThreadPercentage(percent);
    ngcore::SetThreadPercent(percent);
  }

  void BaseStatusHandler::GetStatus (string & str, double & percent)
  {
    // Ng_GetStatus(str, percent);
    ngcore::GetStatus(str, percent);
  }

  void BaseStatusHandler::SetTerminate(void)
  {
    // Ng_SetTerminate();
    // ngcore::SetTerminate();
    ngcore::multithread.terminate = 1;    
  }

  void BaseStatusHandler::UnSetTerminate(void)
  {
    // Ng_UnSetTerminate();
    // ngcore::UnSetTerminate();
    ngcore::multithread.terminate = 1;    
  }

  bool BaseStatusHandler::ShouldTerminate(void)
  {
    // return (Ng_ShouldTerminate() != 0);
    return ngcore::multithread.terminate != 0;
  }

} // namespace ngstd

