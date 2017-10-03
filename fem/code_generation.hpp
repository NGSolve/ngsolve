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
  template <typename T>
  inline string ToLiteral(const T & val)
  {
      stringstream ss;
#if (defined __cpp_hex_float) && (__cpp_hex_float <= __cplusplus)
      ss << std::hexfloat;
      ss << val;
      ss << " /* (" << std::setprecision(16) << std::scientific;
      ss << val << ") */";
#else
      ss << std::setprecision(16);
      ss << val;
#endif
      return ss.str();
  }

  template<>
  inline string ToLiteral (const int &val)
  {
    stringstream ss;
    ss << val;
    return ss.str();
  }


  struct Code
  {
    string top;
    string header;
    string body;
    bool is_simd;
    int deriv;

    string pointer;

    string AddPointer(const void *p );

    static atomic<unsigned> id_counter;
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
    CodeExpr operator ()(int i) { return CodeExpr( S() + '(' + ToLiteral(i) + ')' ); }
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
      return type + " " + code + "("+ToLiteral(value)+");\n";
    }
  };

  inline CodeExpr Var(double val)
  {
    return ToLiteral(val);
  }

  inline CodeExpr Var(Complex val)
  {
    return ToLiteral(val);
  }

  inline CodeExpr Var(string name, int i, int j=0, int k=0)
  {
    return CodeExpr(name + '_' + ToLiteral(i) + '_' + ToLiteral(j) + '_' + ToLiteral(k));
  }

  inline CodeExpr Var(int i, int j=0, int k=0)
  {
    return CodeExpr("var_" + ToLiteral(i) + '_' + ToLiteral(j) + '_' + ToLiteral(k));
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
#else // WIN32
    void *lib;
#endif // WIN32

    void *GetRawFunction( string func_name );

  public:
    Library() : lib(nullptr) {}
    // Compile a given string and load the library
    void Compile( const std::vector<string> &codes );

    void Load( string alib_name );

    template <typename TFunc>
    TFunc GetFunction( string func_name )
    {
      return reinterpret_cast<TFunc>(GetRawFunction(func_name));
    }

    ~Library();

  };
}
#endif // FILE_NGS_CODE_GENERATION___
