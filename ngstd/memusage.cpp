#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#elif defined(__APPLE__)
#include <sys/types.h>
#include <sys/sysctl.h>
#else
#include <unistd.h>
#endif

#include <ngstd.hpp>

namespace ngstd
{
  AllocationException :: AllocationException (const string & s) : Exception(s) { }
  AllocationException :: ~AllocationException () { }

  size_t PhysicalMemory ()
  {
    static size_t mem = [] () -> size_t
    {
#ifdef _WIN32
      MEMORYSTATUSEX status;
      status.dwLength = sizeof(status);
      if (GlobalMemoryStatusEx(&status)) return status.ullTotalPhys;
#elif defined(__APPLE__)
      uint64_t memsize = 0;
      size_t len = sizeof(memsize);
      if (sysctlbyname("hw.memsize", &memsize, &len, nullptr, 0) == 0) return memsize;
#elif defined(_SC_PHYS_PAGES)
      long pages = sysconf(_SC_PHYS_PAGES), pagesize = sysconf(_SC_PAGE_SIZE);
      if (pages > 0 && pagesize > 0) return size_t(pages) * size_t(pagesize);
#endif
      return std::numeric_limits<size_t>::max();
    } ();
    return mem;
  }

  void CheckAllocationSize (size_t n, size_t elsize, const char * what)
  {
    if (elsize && n > std::numeric_limits<size_t>::max() / elsize)
      throw AllocationException (string(what) + ": requested size overflows (" + ToString(n)
                                 + " x " + ToString(elsize) + " bytes)");
    if (n * elsize > PhysicalMemory())
      throw AllocationException (string(what) + ": requested " + ToString(n*elsize)
                                 + " bytes, more than the physical memory ("
                                 + ToString(PhysicalMemory()) + " bytes)");
  }
}
