/*********************************************************************/
/* File:   gpuwrapper.cpp                                            */
/* Author: Joachim Schoeberl                                         */
/* Date:   29. Aug. 2025                                             */
/*********************************************************************/

#include "gpuwrapper.hpp"
#include <stdexcept>
#include <cctype>

namespace ngs_gpu
{
  ArgType ArgType :: FromName (const string & aname)
  {
    // normalize: drop const, collapse whitespace
    string name;
    for (size_t i = 0; i < aname.size(); )
      {
        if (aname.compare(i, 5, "const") == 0 &&
            (i+5 == aname.size() || !std::isalnum((unsigned char)aname[i+5])) &&
            (i == 0 || !std::isalnum((unsigned char)aname[i-1])))
          { i += 5; continue; }
        if (std::isspace((unsigned char)aname[i]))
          {
            if (!name.empty() && name.back() != ' ') name += ' ';
            while (i < aname.size() && std::isspace((unsigned char)aname[i])) i++;
            continue;
          }
        name += aname[i++];
      }
    while (!name.empty() && name.back() == ' ') name.pop_back();
    if (!name.empty() && name.front() == ' ') name.erase(0, 1);
    // "Complex< float >" -> "Complex<float>"
    for (size_t p; (p = name.find(" <")) != string::npos; ) name.erase(p, 1);
    for (size_t p; (p = name.find("< ")) != string::npos; ) name.erase(p+1, 1);
    for (size_t p; (p = name.find(" >")) != string::npos; ) name.erase(p, 1);

    static const std::map<string, ArgType> table =
      {
        { "float", {ArgType::Float, 4} }, { "double", {ArgType::Float, 8} },
        { "half", {ArgType::Float, 2} },
        { "int", {ArgType::Int, 4} }, { "signed int", {ArgType::Int, 4} },
        { "int32_t", {ArgType::Int, 4} }, { "signed", {ArgType::Int, 4} },
        { "unsigned", {ArgType::UInt, 4} }, { "unsigned int", {ArgType::UInt, 4} },
        { "uint", {ArgType::UInt, 4} }, { "uint32_t", {ArgType::UInt, 4} },
        { "long", {ArgType::Int, 8} }, { "long long", {ArgType::Int, 8} },
        { "int64_t", {ArgType::Int, 8} }, { "ptrdiff_t", {ArgType::Int, 8} },
        { "unsigned long", {ArgType::UInt, 8} }, { "unsigned long long", {ArgType::UInt, 8} },
        { "uint64_t", {ArgType::UInt, 8} }, { "size_t", {ArgType::UInt, 8} },
        { "short", {ArgType::Int, 2} }, { "int16_t", {ArgType::Int, 2} },
        { "unsigned short", {ArgType::UInt, 2} }, { "ushort", {ArgType::UInt, 2} },
        { "uint16_t", {ArgType::UInt, 2} },
        { "char", {ArgType::Int, 1} }, { "signed char", {ArgType::Int, 1} },
        { "int8_t", {ArgType::Int, 1} },
        { "unsigned char", {ArgType::UInt, 1} }, { "uint8_t", {ArgType::UInt, 1} },
        { "bool", {ArgType::Bool, 1} },
        { "Complex<float>", {ArgType::Complex, 8} }, { "Complex<double>", {ArgType::Complex, 16} },
      };
    auto it = table.find(name);
    return (it != table.end()) ? it->second : ArgType();
  }

  string ArgType :: ToString() const
  {
    switch (kind)
      {
      case Float: return "float" + std::to_string(8*size);
      case Int:   return "int" + std::to_string(8*size);
      case UInt:  return "uint" + std::to_string(8*size);
      case Bool:  return "bool";
      case Complex: return "complex" + std::to_string(8*size);   // both parts, as numpy
      default:    return size ? std::to_string(size) + " bytes" : "?";
      }
  }


  namespace
  {
    string StripComments (const string & src)
    {
      string out; out.reserve(src.size());
      for (size_t i = 0; i < src.size(); )
        {
          if (src.compare(i, 2, "//") == 0)
            { while (i < src.size() && src[i] != '\n') i++; }
          else if (src.compare(i, 2, "/*") == 0)
            {
              auto e = src.find("*/", i+2);
              i = (e == string::npos) ? src.size() : e+2;
              out += ' ';
            }
          else out += src[i++];
        }
      return out;
    }

    // inside a preprocessor directive (including its continuation lines)?
    bool InDirective (const string & s, size_t pos)
    {
      size_t ls = s.rfind('\n', pos);
      ls = (ls == string::npos) ? 0 : ls+1;
      while (ls >= 2 && s[ls-2] == '\\')
        {
          size_t prev = s.rfind('\n', ls-2);
          ls = (prev == string::npos) ? 0 : prev+1;
        }
      while (ls < s.size() && (s[ls] == ' ' || s[ls] == '\t')) ls++;
      return ls < s.size() && s[ls] == '#';
    }

    string Trim (const string & s)
    {
      size_t a = 0, b = s.size();
      while (a < b && std::isspace((unsigned char)s[a])) a++;
      while (b > a && std::isspace((unsigned char)s[b-1])) b--;
      return s.substr(a, b-a);
    }

    // split at commas not nested in () <> []
    std::vector<string> SplitTopLevel (const string & s)
    {
      std::vector<string> parts;
      int depth = 0; size_t start = 0;
      for (size_t i = 0; i < s.size(); i++)
        {
          char c = s[i];
          if (c == '(' || c == '<' || c == '[') depth++;
          else if (c == ')' || c == '>' || c == ']') depth--;
          else if (c == ',' && depth == 0)
            { parts.push_back (Trim(s.substr(start, i-start))); start = i+1; }
        }
      parts.push_back (Trim(s.substr(start)));
      return parts;
    }

    // "typedef float real;" and "using real = float;"
    std::map<string,string> ParseTypedefs (const string & s)
    {
      std::map<string,string> alias;
      for (size_t pos = s.find("typedef"); pos != string::npos; pos = s.find("typedef", pos+7))
        {
          auto e = s.find(';', pos);
          if (e == string::npos) break;
          auto parts = Trim(s.substr(pos+7, e-pos-7));
          auto sp = parts.find_last_of(" \t\n");
          if (sp != string::npos)
            alias[Trim(parts.substr(sp+1))] = Trim(parts.substr(0, sp));
        }
      for (size_t pos = s.find("using "); pos != string::npos; pos = s.find("using ", pos+6))
        {
          auto e = s.find(';', pos);
          auto eq = s.find('=', pos);
          if (e == string::npos) break;
          if (eq != string::npos && eq < e)
            alias[Trim(s.substr(pos+6, eq-pos-6))] = Trim(s.substr(eq+1, e-eq-1));
        }
      return alias;
    }

    ArgType ResolveType (string name, const std::map<string,string> & alias)
    {
      for (int i = 0; i < 8; i++)   // alias chains
        {
          auto t = ArgType::FromName (name);
          if (t.Known()) return t;
          auto it = alias.find (Trim(name));
          if (it == alias.end()) return t;
          name = it->second;
        }
      return {};
    }
  }


  std::map<string, KernelSignature> ParseKernelSignatures (const string & source)
  {
    std::map<string, KernelSignature> sigs;
    string src = StripComments (source);
    auto alias = ParseTypedefs (src);

    for (size_t pos = src.find("KERNEL"); pos != string::npos; pos = src.find("KERNEL", pos+6))
      {
        if (pos > 0 && (std::isalnum((unsigned char)src[pos-1]) || src[pos-1] == '_')) continue;
        size_t p = pos+6;
        bool bounds = src.compare(p, 7, "_BOUNDS") == 0;
        if (bounds) p += 7;
        while (p < src.size() && std::isspace((unsigned char)src[p])) p++;
        if (p >= src.size() || src[p] != '(') continue;
        if (InDirective (src, pos)) continue;

        // balanced argument list
        int depth = 0; size_t end = p;
        for (; end < src.size(); end++)
          {
            if (src[end] == '(') depth++;
            else if (src[end] == ')' && --depth == 0) break;
          }
        if (end >= src.size()) break;

        auto parts = SplitTopLevel (src.substr(p+1, end-p-1));
        size_t first = bounds ? 3 : 1;
        if (parts.size() < first) continue;
        string name = parts[0];

        KernelSignature sig;
        for (size_t i = first; i < parts.size(); i++)
          {
            KernelSignature::Arg arg { KernelSignature::Unknown, {}, "" };
            const string & a = parts[i];
            auto lp = a.find('(');
            if (lp != string::npos && a.back() == ')')
              {
                string macro = Trim(a.substr(0, lp));
                auto inner = SplitTopLevel (a.substr(lp+1, a.size()-lp-2));
                if (inner.size() == 2)
                  {
                    string tname = inner[0];
                    if (macro == "GLOBAL" || macro == "GLOBAL_IN" || macro == "GLOBAL_ATOMIC")
                      arg.kind = KernelSignature::Buffer;
                    else if (macro == "GLOBAL_ATOMIC_COMPLEX")
                      {
                        // declared with the real type, holds Complex<T>
                        arg.kind = KernelSignature::Buffer;
                        tname = "Complex<" + Trim(inner[0]) + ">";
                      }
                    else if (macro == "VALUE")
                      arg.kind = KernelSignature::Value;
                    if (arg.kind != KernelSignature::Unknown)
                      {
                        arg.type = ResolveType (tname, alias);
                        arg.name = inner[1];
                      }
                  }
              }
            sig.args.push_back (arg);
          }
        sigs[name] = std::move(sig);
      }
    return sigs;
  }


  void CheckKernelArgs (const Kernel & kernel, const std::vector<KernelArg> & args)
  {
    auto sig = kernel.Signature();
    if (!sig) return;

    auto fail = [&] (const string & msg)
    { throw std::runtime_error ("kernel " + kernel.Name() + ": " + msg); };

    if (sig->args.size() != args.size())
      fail ("expects " + std::to_string(sig->args.size()) + " arguments, got " +
            std::to_string(args.size()));

    for (size_t i = 0; i < args.size(); i++)
      {
        const auto & decl = sig->args[i];
        const auto & arg = args[i];
        if (decl.kind == KernelSignature::Unknown) continue;

        string where = "argument " + std::to_string(i) + " (" + decl.name + ")";
        bool isbuf = arg.GetKind() == KernelArg::Kind::Buffer;
        if (isbuf != (decl.kind == KernelSignature::Buffer))
          fail (where + ": expected " + (isbuf ? "a value, got a buffer" : "a buffer, got a value"));

        ArgType want = decl.type, got = arg.GetType();
        if (!want.Known()) continue;
        string what = isbuf ? "buffer of " : "value of type ";
        if (got.Known() && got != want)
          fail (where + ": expected " + what + want.ToString() + ", got " + got.ToString());
        if (!isbuf && got.size && got.size != want.size)
          fail (where + ": expected " + what + want.ToString() + ", got " + got.ToString());
      }
  }


  static DeviceCreator device_creator;
  static shared_ptr<Device> the_device;

  void SetDeviceCreator (DeviceCreator creator)
  {
    device_creator = std::move(creator);
    the_device.reset();
  }

  bool HasDevice() { return bool(device_creator); }

  shared_ptr<Device> GetDevice()
  {
    if (!the_device && device_creator)
      the_device = device_creator();
    return the_device;
  }
}
