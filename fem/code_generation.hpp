#ifndef FILE_NGS_CODE_GENERATION___
#define FILE_NGS_CODE_GENERATION___

/*********************************************************************/
/* File:   code_generation.hpp                                       */
/* Author: Matthias Hochsteger                                       */
/* Date:   13. Apr. 2016                                             */
/*********************************************************************/

#ifdef WIN32
#include <windows.h>
#else // WIN32
#include <dlfcn.h>
#endif //WIN32

#include <map>

namespace ngfem
{
  struct Code
  {
    string header;
    string body;
    bool is_simd;
    int deriv;

    static string Map( string code, std::map<string,string> variables ) {
      for ( auto mapping : variables ) {
        string oldStr = '{'+mapping.first+'}';
        string newStr = mapping.second;
        size_t pos = 0;
        while((pos = code.find(oldStr, pos)) != std::string::npos){
          code.replace(pos, oldStr.length(), newStr);
          pos += newStr.length();

        }
      }
      return code;
    }
  };

  struct CodeExpr
  {
    string code;
    CodeExpr( string acode="" ) : code(acode) {;}
    string S() const { return code; }

    string Op(char c) { return code.size() ? string(" ") + c + ' ' : ""; }

    CodeExpr operator +(CodeExpr other) { return CodeExpr(string("(") + S()+Op('+')+other.S() + ')'); }
    CodeExpr operator -(CodeExpr other) { return CodeExpr(string("(") + S()+Op('-')+other.S() + ')'); }
    CodeExpr operator *(CodeExpr other) { return CodeExpr(string("(") + S()+Op('*')+other.S() + ')'); }
    CodeExpr operator /(CodeExpr other) { return CodeExpr(string("(") + S()+Op('/')+other.S() + ')'); }
    void operator +=(CodeExpr other) { code = "(" + S()+Op('+')+other.S() + ')'; }
    void operator -=(CodeExpr other) { code = "(" + S()+Op('-')+other.S() + ')'; }
    void operator *=(CodeExpr other) { code = "(" + S()+Op('*')+other.S() + ')'; }
    void operator /=(CodeExpr other) { code = "(" + S()+Op('/')+other.S() + ')'; }

    operator string () { return code; }
    CodeExpr operator ()(int i) { return CodeExpr( S() + '(' + ToString(i) + ')' ); }
    CodeExpr Func(string s) { return CodeExpr( s + "(" + S() + ")" ); }
    CodeExpr Call(string s, string args="") { return CodeExpr( S()+'.'+ s + "(" + args + ")"); }
    string Assign (CodeExpr other, bool declare = true)
    {
      string result;
      if(declare)
        result += "auto ";
      result += S()+" = "+other.S() + ";\n";
      return result;
    }

    string Declare(string type = "auto")
    {
      return type + " " + code + ";\n";
    }

    template<typename TVal>
    string Declare(string type, TVal value )
    {
      return type + " " + code + "("+ToString(value)+");\n";
    }
  };

  inline CodeExpr Var(double val)
  {
    return CodeExpr(ToString(val));
  }

  inline CodeExpr Var(Complex val)
  {
    string res("Complex(");
    res += ToString(val.real());
    res += ",";
    res += ToString(val.imag());
    res += ")";
    return res;
  }

  inline CodeExpr Var(string name, int i, int j=0, int k=0)
  {
    return CodeExpr(name + '_' + ToString(i) + '_' + ToString(j) + '_' + ToString(k));
  }

  inline CodeExpr Var(int i, int j=0, int k=0)
  {
    return CodeExpr("var_" + ToString(i) + '_' + ToString(j) + '_' + ToString(k));
  }

  template<typename TFunc>
  void TraverseDimensions( FlatArray<int> dims, const TFunc &func)

  {
    switch(dims.Size())
    {
      case 0:
        func(0,0,0);
        break;
      case 1:
        for (int i : Range(max2(1, dims[0])))
          func(i,i,0);
        break;
      case 2:
        for (int i : Range(max2(1, dims[0])))
          for (int j : Range(max2(1, dims[1])))
            func(i*dims[1]+j, i, j);
        break;
      default:
        throw Exception("TraverseDimensions: too many dimensions!");
    }
  }

  static void GetIndex( FlatArray<int> dims, int i, int &iout, int &jout )
  {
    iout = jout = 0;
    switch(dims.Size())
    {
      case 0:
        break;
      case 1:
        iout = i;
        break;
      case 2:
        iout = i/dims[1];
        jout = i%dims[1];
        break;
      default:
        throw Exception("GetIndex: too many dimensions!");
    }
  }

  class Library
  {
    static int counter; // defined in coefficient.cpp
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
    Library() : lib(nullptr) {}
    // Compile a given string and load the library
    void Compile( string code )
    {
      static ngstd::Timer tcompile("CompiledCF::Compile");
      static ngstd::Timer tlink("CompiledCF::Link");
      string file_prefix = "code" + ToString(counter++);
      ofstream codefile(file_prefix+".cpp");
      codefile << code;
      codefile.close();
      cout << IM(3) << "compiling..." << endl;
      tcompile.Start();
#ifdef WIN32
      string scompile = "cmd /C \"ngscxx.bat " + file_prefix + ".cpp\"";
      string slink = "cmd /C \"ngsld.bat " + file_prefix + ".obj\"";
#else
      string scompile = "ngscxx -c " + file_prefix + ".cpp -o " + file_prefix + ".o";
      string slink = "ngsld -shared " + file_prefix + ".o -o " + file_prefix + ".so -lngstd -lngfem";
#endif
      int err = system(scompile.c_str());
      if (err) throw Exception ("problem calling compiler");
      tcompile.Stop();
      cout << IM(3) << "linking..." << endl;
      tlink.Start();
      err = system(slink.c_str());
      if (err) throw Exception ("problem calling linker");      
      tlink.Stop();
      cout << IM(3) << "done" << endl;
#ifdef WIN32
      Load(file_prefix+".dll");
#else
      Load(file_prefix+".so");
#endif
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
        if(rc != 0) cerr << "Failed to close library " << lib_name << endl;
      }
    }
#endif // WIN32
  };
}
#endif // FILE_NGS_CODE_GENERATION___
