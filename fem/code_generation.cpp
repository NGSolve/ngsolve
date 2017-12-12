#include "fem.hpp"
#include <algorithm>

namespace ngfem
{
    atomic<unsigned> Code::id_counter{0};

    void Code::AddLinkFlag(string flag)
    {
        if(std::find(std::begin(link_flags), std::end(link_flags), flag) == std::end(link_flags))
          link_flags.push_back(flag);
    }

    string Code::AddPointer(const void *p)
    {
        string name = "compiled_code_pointer" + ToString(id_counter++); 
        top += "extern \"C\" void* " + name + ";\n";
        stringstream s_ptr;
#ifdef WIN32
        s_ptr << "0x";
#endif
        s_ptr << std::hex << p;
        pointer += "void *" + name + " = reinterpret_cast<void*>(" + s_ptr.str() + ");\n";
        return name;
    }

    unique_ptr<SharedLibrary> CompileCode(const std::vector<string> &codes, const std::vector<string> &link_flags )
    {
      static int counter = 0;
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
        string scompile = "ngscxx -c " + file_prefix + ".cpp -o " + file_prefix + ".o";
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
        string slink = "ngsld -shared " + object_files + " -o " + prefix + ".so -lngstd -lngbla -lngfem";
        for (auto flag : link_flags)
            slink += " "+flag;
#endif
      int err = system(slink.c_str());
      if (err) throw Exception ("problem calling linker");      
      tlink.Stop();
      cout << IM(3) << "done" << endl;
      auto library = make_unique<SharedLibrary>();
#ifdef WIN32
      library->Load(prefix+".dll");
#else
      library->Load("./"+prefix+".so");
#endif
      return library;
    }

}
