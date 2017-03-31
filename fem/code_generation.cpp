#include "fem.hpp"

namespace ngfem
{
    unsigned Code::id_counter = 0;

    string Code::AddPointer(const void *p)
    {
        string name = "compiled_code_pointer" + ToString(id_counter++); 
        top += "extern \"C\" void* " + name + ";\n";
        stringstream s_ptr;
        s_ptr << std::hex << p;
        pointer += "void *" + name + " = reinterpret_cast<void*>(0x" + s_ptr.str() + ");\n";
        return name;
    }

    void Library::Load( string alib_name )
    {
      lib_name = alib_name;
#ifdef WIN32
      lib = LoadLibrary(lib_name.c_str());
      if (!lib) throw std::runtime_error(string("Could not load library ") + lib_name);
#else // WIN32
      lib = dlopen(lib_name.c_str(), RTLD_NOW);
      if(lib == nullptr) throw std::runtime_error(dlerror());
#endif // WIN32
    }

    void Library::Compile( std::vector<string> &codes )
    {
      static ngstd::Timer tcompile("CompiledCF::Compile");
      static ngstd::Timer tlink("CompiledCF::Link");
      string object_files;
      int i = 0;
      string prefix = "code" + ToString(counter++);
      for(string code : codes) {
        string file_prefix = prefix+"_"+ToString(i++);
        ofstream codefile(file_prefix+".cpp");
        codefile << code;
        codefile.close();
        cout << IM(3) << "compiling..." << endl;
        tcompile.Start();
#ifdef WIN32
        string scompile = "cmd /C \"ngscxx.bat " + file_prefix + ".cpp\"";
        object_files += file_prefix+".obj ";
#else
        string scompile = "ccache ngscxx -c " + file_prefix + ".cpp -o " + file_prefix + ".o";
        object_files += file_prefix+".o ";
#endif
        int err = system(scompile.c_str());
        if (err) throw Exception ("problem calling compiler");
        tcompile.Stop();
      }

      cout << IM(3) << "linking..." << endl;
      tlink.Start();
#ifdef WIN32
        string slink = "cmd /C \"ngsld.bat /OUT:" + prefix+".dll " + object_files + "\"";
#else
        string slink = "ngsld -shared " + object_files + " -o " + prefix + ".so -lngstd -lngfem";
#endif
      int err = system(slink.c_str());
      if (err) throw Exception ("problem calling linker");      
      tlink.Stop();
      cout << IM(3) << "done" << endl;
#ifdef WIN32
      Load(prefix+".dll");
#else
      Load(prefix+".so");
#endif
    }

    void* Library::GetRawFunction( string func_name )
    {
#ifdef WIN32
      void* func = GetProcAddress(lib, func_name.c_str());
      if(func == nullptr)
        throw std::runtime_error(string("Could not find function ") + func_name + " in library " + lib_name);

#else // WIN32
      void* func = dlsym(lib, func_name.c_str());
      if(func == nullptr)
          throw std::runtime_error(dlerror());

#endif // WIN32

      return func;
    }

    Library::~Library()
    {
      if(lib)
      {

#ifdef WIN32
        FreeLibrary(lib);
#else // WIN32
        int rc = dlclose(lib);
        if(rc != 0) cerr << "Failed to close library " << lib_name << endl;
#endif // WIN32
      }
    }
}
