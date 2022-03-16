#include "fem.hpp"
#include <algorithm>
#include<l2hofe_impl.hpp>
#include<l2hofefo.hpp>
#include<regex>
#include<cstdio>

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

    unique_ptr<SharedLibrary> CompileCode(const std::vector<std::variant<filesystem::path, string>> &codes, const std::vector<string> &link_flags, bool keep_files )
    {
      static int counter = 0;
      static ngstd::Timer tcompile("CompiledCF::Compile");
      static ngstd::Timer tlink("CompiledCF::Link");
      string object_files;
      int rank = 0;
#ifdef PARALLEL
      rank = ngcore::NgMPI_Comm(MPI_COMM_WORLD).Rank();
#endif // PARALLEL
      filesystem::path lib_dir;
#ifdef WIN32
      lib_dir = filesystem::path(std::tmpnam(nullptr)).concat("_ngsolve_"+ToString(rank)+"_"+ToString(counter++));
#else // WIN32
      string tmp_template = filesystem::temp_directory_path().append("ngsolve_tmp_"+ToString(rank)+"_"+ToString(counter++)+"_XXXXXX");
      if(mkdtemp(&tmp_template[0])==nullptr)
          throw Exception("could not create temporary directory");

      lib_dir = tmp_template;
#endif // WIN32
      string chdir_cmd = "cd " + lib_dir.string() + " && ";
      filesystem::create_directories(lib_dir);
      for(auto i : Range(codes.size())) {
        filesystem::path src_file;
        if(std::holds_alternative<filesystem::path>(codes[i]))
            src_file = filesystem::absolute(std::get<filesystem::path>(codes[i]));
        else
        {
            string code = std::get<string>(codes[i]);
            src_file = filesystem::path(lib_dir).append("code_" + ToString(i) + ".cpp");
            ofstream codefile(src_file);
            codefile << code;
            codefile.close();
        }
        cout << IM(3) << "compiling..." << endl;
        tcompile.Start();
#ifdef WIN32
        string scompile = "cmd /C \"ngscxx.bat " + src_file.string();
        object_files += " " + filesystem::path(src_file).replace_extension(".obj ").string();
#else // WIN32
        auto obj_file = " " + filesystem::path(src_file).replace_extension(".o").string();
        string scompile = "ngscxx -c " + src_file.string() + " -o " + obj_file;
        object_files += obj_file;
#endif // WIN32
        int err = system((chdir_cmd + scompile).c_str());
        if (err) throw Exception ("problem calling compiler");
        tcompile.Stop();
      }

      cout << IM(3) << "linking..." << endl;
      tlink.Start();
      auto lib_file = filesystem::path(lib_dir).append("library");
#ifdef WIN32
      lib_file.concat(".dll");
      string slink = "cmd /C \"ngsld.bat /OUT:" + lib_file.string() + " " + object_files;
      for (auto flag : link_flags)
        slink += " "+flag;
      slink += " \"";
#else // WIN32
      lib_file.concat(".so");
      string slink = "ngsld -shared " + object_files + " -o " + lib_file.string() + " -lngstd -lngbla -lngfem -lngcomp -lngcore";
      for (auto flag : link_flags)
        slink += " "+flag;
#endif // WIN32
      int err = system((chdir_cmd + slink).c_str());
      if (err) throw Exception ("problem calling linker");      
      tlink.Stop();
      cout << IM(3) << "done" << endl;
      if(keep_files)
      {
          cout << IM(2) << "keeping generated files at " << lib_dir.string() << endl;
          return make_unique<SharedLibrary>(lib_file);
      }
      else
          return make_unique<SharedLibrary>(lib_file, lib_dir);
    }

    namespace detail {
        // T_CalcShape is protected, thus we have to derive
        template <ELEMENT_TYPE ET, typename BASE>
          class MyFEL : public BASE {
            public:
              using BASE::ndof;
              using BASE::vnums;
              using BASE::order;
              using BASE::order_inner;
              using BASE::GetFaceSort;
              using BASE::GetEdgeSort;
              MyFEL(int order) : BASE(order) {}
              MyFEL() : BASE() {}
              template<typename Tx, typename TFA>
                INLINE void MyCalcShape (TIP<ET_trait<ET>::DIM,Tx> ip, TFA & shape) const
                {
                    BASE::T_CalcShape(ip, shape);
                }
              void GetDiagMassMatrix(FlatVector<> mass) const {
                  throw Exception("GetDiagMassMatrix not implemented");
              }
          };

        // This class behaves like a numeric type, but 'records' all computations in expression strings
        struct CCode {
            static std::vector<string> expressions;
            static int find( std::vector<string> &v, string val ){
                int i = 0;
                for(auto &s : v) {
                    if(s==val)
                        return i;
                    i++;
                }
                return -1;
            }

            static string strip(string s) {
                int n = s.size();
                if(n<=1) return s;
                if(s[0] == '(' && s[n-1] == ')') return strip(s.substr(1,n-2));
                return s;
            }
            mutable string s;

            static CCode Par(const CCode &c) {
                return c.s;
            }

            void Check() {
                static string int_num = "var[0-9]*";
                static regex pattern(int_num);
                if(s=="") return;
                int index = find( expressions, strip(s));
                if(index>=0) {
                    s = "var"+ToString(index);
                }
                else {
                    if(!regex_match(strip(s), pattern)) {
                        expressions.push_back(strip(s));
                        s = "var"+ToString(expressions.size()-1);
                    }
                }
            }

            CCode(const CCode &c) :
              s(c.s)
            {
                Check();
            }

            CCode(string as = "") : s(as) {
                Check();
            }

            CCode(double val) {
                std::stringstream str;
                str << fixed << setprecision( 15 ) << val;
                s = str.str();
                Check();
            }

            virtual CCode operator +(const CCode &c) { return CCode(s+'+'+c.s); }
            virtual CCode operator -(const CCode &c) { return CCode(s+'-'+c.s); }
            virtual CCode operator -() { return CCode('-'+s); }
            virtual CCode operator *(const CCode &c) { return CCode(s+'*'+c.s); }
            virtual void operator +=(const CCode &c) { *this = *this+c; }
            virtual void operator *=(const CCode &c) { *this = *this*c; }
            virtual CCode &operator=(const CCode &c) {
                s = c.s;
                return *this;
            }
            virtual CCode operator /(const CCode &c) { return CCode(s+'/'+c.s); }
        };

        std::vector<string> CCode::expressions;
        CCode operator -(int val, const CCode &c) { return CCode(1.0*val)-c; }
        CCode operator *(int val, const CCode &c) { return CCode(1.0*val)*c; }
        CCode operator -(double val, const CCode &c) { return CCode(val)-c; }
        CCode operator *(double val, const CCode &c) { return CCode(val)*c; }

        ostream &operator <<(ostream & s, const CCode &c) {
            s << c.s;
            return s;
        }

        // Generate code for L2HighOrder elements
        string GenerateL2ElementCode(int order) {
            auto genCode = [&] (const auto &fel, auto ip, string elname) -> string
            {
                stringstream f;
                f <<
                  "float Eval" << elname << "(int element, float x, float y, float z )\n"
                  "{                             \n"
                  " float result = 0.0;" << endl;

                stringstream ss;
                fel.MyCalcShape (ip, SBLambda([&] (int i, auto c) {
                                              ss << "result += texelFetch( coefficients, element*"+ToString(fel.ndof) + "+"  + ToString(i) + ").r * " + c.s << ";" << endl;
                                              }));

                int i = 0;
                for(auto &s : CCode::expressions)
                    f << "float var" << ToString(i++) << " = " <<  s << ";" << endl;
                f << ss.str() << endl;
                f << "return result;" << endl;
                f << "}" << endl;
                return f.str();
            };

            string code;
            {
              CCode::expressions.clear(); CCode x("x"); CCode y("y"); TIP<2,CCode> ip(x,y, -1, VOL);
                code += genCode(MyFEL<ET_TRIG, L2HighOrderFE<ET_TRIG>>(order), ip, "TRIG");
            } {
              CCode::expressions.clear(); CCode x("x"); CCode y("y"); TIP<2,CCode> ip(x,y, -1, VOL);
                code += genCode(MyFEL<ET_QUAD, L2HighOrderFE<ET_QUAD>>(order), ip, "QUAD");
            } {
              CCode::expressions.clear(); CCode x("x"); CCode y("y"); CCode z("z"); TIP<3,CCode> ip(x,y,z, -1, VOL);
                code += genCode(MyFEL<ET_TET, L2HighOrderFE<ET_TET>>(order), ip, "TET");
            } {
              CCode::expressions.clear(); CCode x("x"); CCode y("y"); CCode z("z"); TIP<3,CCode> ip(x,y,z, -1, VOL);
                code += genCode(MyFEL<ET_HEX, L2HighOrderFE<ET_HEX>>(order), ip, "HEX");
            } {
              CCode::expressions.clear(); CCode x("x"); CCode y("y"); CCode z("z"); TIP<3,CCode> ip(x,y,z, -1, VOL);
                code += genCode(MyFEL<ET_PYRAMID, L2HighOrderFE<ET_PYRAMID>>(order), ip, "PYRAMID");
            } {
              CCode::expressions.clear(); CCode x("x"); CCode y("y"); CCode z("z"); TIP<3,CCode> ip(x,y,z, -1, VOL);
                code += genCode(MyFEL<ET_PRISM, L2HighOrderFE<ET_PRISM>>(order), ip, "PRISM");
            }
            return code;
        }

    }
}
