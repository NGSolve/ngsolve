#ifdef WIN32
#include <windows.h>
#else // WIN32
#include <dlfcn.h>
#endif //WIN32

#include <string>

using std::string;

#ifdef WIN32
class Library
{
  HINSTANCE lib;
  string lib_name;
public:
  void Load( string alib_name )
    {
      lib_name = alib_name;
      lib = LoadLibrary(lib_name.c_str());
      if (!lib) throw std::runtime_error(string("Could not load library ") + lib_name);
    }

  template <typename TFunc>
    TFunc GetFunction( string func_name )
      {
        TFunc func = reinterpret_cast<TFunc>(GetProcAddress(lib, func_name.c_str()));
        if(func == nullptr)
          throw std::runtime_error(string("Could not find function ") + func_name + " in library " + lib_name);
        return func;
      }

  ~Library()
    {
      if(lib)
        FreeLibrary(lib);
    }
};

#else // WIN32
class Library
{
  void *lib;
  string lib_name;

public:
  void Load( string alib_name )
    {
      lib_name = alib_name;
      lib = dlopen(lib_name.c_str(), RTLD_NOW);
      if(lib == nullptr) throw std::runtime_error(dlerror());
    }

  template <typename TFunc>
    TFunc GetFunction( string func_name )
      {

        TFunc func = reinterpret_cast<TFunc>(dlsym(lib, func_name.c_str()));
        if(func == nullptr) throw std::runtime_error(dlerror());
        return func;
      }

  ~Library()
    {
      if(lib)
        {
          int rc = dlclose(lib);
          if(rc != 0) throw std::runtime_error(string("Failed to close ") + lib_name);
        }
    }
};
#endif // WIN32
