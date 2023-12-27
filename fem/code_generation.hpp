#ifndef FILE_NGS_CODE_GENERATION___
#define FILE_NGS_CODE_GENERATION___

/*********************************************************************/
/* File:   code_generation.hpp                                       */
/* Author: Matthias Hochsteger                                       */
/* Date:   13. Apr. 2016                                             */
/*********************************************************************/

#include <filesystem>
#include <map>
#include <variant>

#include <core/ngcore.hpp>

namespace ngfem
{
  using namespace ngbla;
  
  NGS_DLL_HEADER extern bool code_uses_tensors;

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
      ss << std::setprecision(16) << std::scientific;
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

    string res_type;
    bool is_simd;
    int deriv;
    std::vector<string> link_flags;

    string pointer;

    NGS_DLL_HEADER string AddPointer(const void *p );

    void AddLinkFlag(string flag);

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

    [[deprecated("use Declare(i,dims,iscomplex) instead")]]  
    void Declare (string type, int i, FlatArray<int> dims);
    void Declare (int i, FlatArray<int> dims, bool iscomplex);
    string GetType (bool iscomplex) const;
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
    CodeExpr operator ()(int i, int j) { return CodeExpr( S() + '(' + ToLiteral(i) + ',' + ToLiteral(j) + ')' ); }
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
    return "Complex"+ToLiteral(val);
  }

  /*
  inline CodeExpr Var(string name, int i, int j=0, int k=0)
  {
    return CodeExpr(name + '_' + ToLiteral(i) + '_' + ToLiteral(j) + '_' + ToLiteral(k));
  }
  */
  inline CodeExpr Var(string name, int i)
  {
    /*
    if (code_uses_tensors)
      return CodeExpr(name + '_' + ToLiteral(i) + "()");
    else
    */
      return CodeExpr(name + '_' + ToLiteral(i));
  }
  inline CodeExpr Var(string name, int i, int j)
  {
    /*
    if (code_uses_tensors)
      return CodeExpr(name + '_' + ToLiteral(i) + '(' + ToLiteral(j) + ')');
    else
    */
      return CodeExpr(name + '_' + ToLiteral(i) + '_' + ToLiteral(j));
  }
  inline CodeExpr Var(string name, int i, int j, int k)
  {
    /*
    if (code_uses_tensors)
      return CodeExpr(name + '_' + ToLiteral(i) + '(' + ToLiteral(j) + ',' + ToLiteral(k) + ')');
    else
    */
      return CodeExpr(name + '_' + ToLiteral(i) + '_' + ToLiteral(j) + '_' + ToLiteral(k));
  }

  // linear index of tensor of dimensions dims
  inline CodeExpr Var(string name, int i, int index, FlatArray<int> dims)
  {
    ArrayMem<int,8> ind(dims.Size());
    for (int j = dims.Size()-1; j >= 0; j--)
      {
        ind[j] = index % dims[j];
        index /= dims[j];
      }
    /*
    if (code_uses_tensors)
      {
        string str = name + '_' + ToLiteral(i) + "(";
        for (int j = 0; j < ind.Size(); j++)
          {
            if (j > 0) str += ',';
            str += ToLiteral(ind[j]);
          }
        str += ")";
        return CodeExpr(str);
      }
    else
    */
      {
        string str = name + '_' + ToLiteral(i);
        for (int j = 0; j < ind.Size(); j++)
          str += '_' + ToLiteral(ind[j]);
        return CodeExpr(str);
      }
  }


  
  /*
  inline CodeExpr Var(int i, int j=0, int k=0)
  {
    return CodeExpr("var_" + ToLiteral(i) + '_' + ToLiteral(j) + '_' + ToLiteral(k));
  }
  */
  inline CodeExpr Var(int i)
  {
    if (code_uses_tensors)
      return CodeExpr("var_" + ToLiteral(i) + "()");
    else
      return CodeExpr("var_" + ToLiteral(i));
  }

  inline CodeExpr Var(int i, int j)
  {
    if (code_uses_tensors)
      return CodeExpr("var_" + ToLiteral(i) + '(' + ToLiteral(j) + ')');
    else
      return CodeExpr("var_" + ToLiteral(i) + '_' + ToLiteral(j));
  }

  inline CodeExpr Var(int i, int j, int k)
  {
    if (code_uses_tensors)
      return CodeExpr("var_" + ToLiteral(i) + '(' + ToLiteral(j) + ',' + ToLiteral(k) + ')');
    else
      return CodeExpr("var_" + ToLiteral(i) + '_' + ToLiteral(j) + '_' + ToLiteral(k));
  }
  
  // linear index of tensor of dimensions dims
  inline CodeExpr Var(int i, int index, FlatArray<int> dims)
  {
    ArrayMem<int,8> ind(dims.Size());
    for (int j = dims.Size()-1; j >= 0; j--)
      {
        ind[j] = index % dims[j];
        index /= dims[j];
      }

    if (code_uses_tensors)
      {
        string str = "var_" + ToLiteral(i) + "(";
        for (int j = 0; j < ind.Size(); j++)
          {
            if (j > 0) str += ',';
            str += ToLiteral(ind[j]);
          }
        str += ")";
        return CodeExpr(str);
      }
    else
      {
        string str = "var_" + ToLiteral(i);
        for (int j = 0; j < ind.Size(); j++)
          str += '_' + ToLiteral(ind[j]);
        return CodeExpr(str);
      }
      // return CodeExpr("var_" + ToLiteral(i) + '_' + ToLiteral(j) + '_' + ToLiteral(k));
  }

  
  inline void GetIndex( FlatArray<int> dims, int i, int &iout, int &jout )
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
        {
          int d1 = dims[1];
          iout = i/d1;
          jout = i%d1;
          break;
        }
      default:
        throw Exception("GetIndex: too many dimensions!");
    }
  }

  std::filesystem::path CreateTempDir();
  unique_ptr<SharedLibrary> CompileCode(const std::vector<std::variant<filesystem::path, string>> &codes, const std::vector<string> &link_flags, bool keep_files = false );
  namespace detail {
      string GenerateL2ElementCode(int order);
  }
  using detail::GenerateL2ElementCode;
}
#endif // FILE_NGS_CODE_GENERATION___
