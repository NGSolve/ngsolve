#include "fem.hpp"
#include <algorithm>
#include<l2hofe_impl.hpp>
#include<l2hofefo.hpp>
#include<regex>

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
        string slink = "ngsld -shared " + object_files + " -o " + prefix + ".so -lngstd -lngbla -lngfem -lngcomp -lngcore";
#endif
      for (auto flag : link_flags)
        slink += " "+flag;
      int err = system(slink.c_str());
      if (err) throw Exception ("problem calling linker");      
      tlink.Stop();
      cout << IM(3) << "done" << endl;
      auto library = make_unique<SharedLibrary>();
#ifdef WIN32
      library->Load(prefix+".dll");
#else
      char *temp = getcwd(nullptr, 0);
      string cwd(temp);
      free(temp);
      library->Load(cwd+"/"+prefix+".so");
#endif
      return library;
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
