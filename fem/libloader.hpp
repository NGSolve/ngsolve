#ifdef WIN32
#include <windows.h>
#else // WIN32
#include <dlfcn.h>
#endif //WIN32

#include <string>
#include <ngstd.hpp>

using std::string;

class Library
{
  static int counter;
  string lib_name;

#ifdef WIN32
  HINSTANCE lib;
  void Load( string alib_name )
    {
      lib_name = alib_name;
      lib = LoadLibrary(lib_name.c_str());
      if (!lib) throw std::runtime_error(string("Could not load library ") + lib_name);
    }
#else // WIN32
  void *lib;
  void Load( string alib_name )
    {
      lib_name = alib_name;
      lib = dlopen(lib_name.c_str(), RTLD_NOW);
      if(lib == nullptr) throw std::runtime_error(dlerror());
    }
#endif // WIN32


public:
  // Compile a given string and load the library
  void Compile( string code ) {
      static ngstd::Timer tcompile("CompiledCF::Compile");
      static ngstd::Timer tlink("CompiledCF::Link");
      string file_prefix = "code" + ToString(counter);
      ofstream codefile(file_prefix+".cpp");
      codefile << code;
      codefile.close();
//       string edit_str = "vim " + file_prefix+".cpp";
//       system(edit_str.c_str());
      cout << "compiling..." << endl;
      tcompile.Start();
      string scompile = "ngscxx -c " + file_prefix + ".cpp -o " + file_prefix + ".o";
      system(scompile.c_str());
      tcompile.Stop();
//       scompile = "ngscxx -S -g0 " + file_prefix + ".cpp -o " + file_prefix + ".s";
//       system(scompile.c_str());
//       scompile = "./cl -S -g0 " + file_prefix + ".cpp -o " + file_prefix + "_clang.s";
//       system(scompile.c_str());
      cout << "linking..." << endl;
      tlink.Start();
      string slink = "ngsld -shared " + file_prefix + ".o -o " + file_prefix + ".so -lngstd -lngfem";
      system(slink.c_str());
      tlink.Stop();
      cout << "done" << endl;
      Load(file_prefix+".so");
      counter++;
  }


#ifdef WIN32
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
#else // WIN32
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
#endif // WIN32
};
