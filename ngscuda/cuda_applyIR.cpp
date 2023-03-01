
#include <la.hpp>
#include <comp.hpp>
#include "cuda_linalg.hpp"

using namespace ngcomp;

namespace ngla
{
  
  class DevApplyIntegrationPoints : public DevMatrix
  {
    size_t h, w;
    size_t dimx, dimy, nip;
    unique_ptr<SharedLibrary> library;

    typedef void (*lib_function)(size_t nip, double * input, size_t dist_input,
                                 double * output, size_t dist_output);
    lib_function compiled_function = nullptr;
    
  public:
    DevApplyIntegrationPoints (const ApplyIntegrationPoints & aipmat)
    {
      h = aipmat.Height();
      w = aipmat.Width();

      dimx = aipmat.GetDimX();
      dimy = aipmat.GetDimY();
      nip = aipmat.GetNIP();
      
      // generate cuda code, similar as for host (with C-Function + Kernel):

      auto & trialproxies = aipmat.GetTrialProxies();
      
      Array<int> proxyoffset;
      int starti = 0;
      for (auto proxy : trialproxies)
        {
          proxyoffset.Append (starti);
          starti += proxy->Evaluator()->Dim();
        }
      
      stringstream s;
      s <<
        "#include <cstddef>\n"
        "__global__ void ApplyIPFunctionKernel (size_t nip, double * input, size_t dist_input,\n"
        "                      double * output, size_t dist_output) {\n";
      
      int base_output = 0;
      for (auto cf : aipmat.GetCFs())
        {
          auto compiledcf = Compile (cf, false);
          Code code = compiledcf->GenerateProgram(0, false);
          
          s << "{\n";
          // cout << code.header << endl;
          
          for (auto step : Range(compiledcf->Steps()))
            if (auto proxycf = dynamic_cast<ProxyFunction*> (compiledcf->Steps()[step]))
              if (auto pos = trialproxies.Pos(proxycf); pos != trialproxies.ILLEGAL_POSITION)
                {
                  s << "auto values_" << step << " = [dist_input,input](size_t i, int comp)\n"
                    " { return input[i + (comp+" << proxyoffset[pos] << ")*dist_input]; };\n";
                  s << "bool constexpr has_values_" << step << " = true;\n" << endl;
                }


          s << "int tid = blockIdx.x*blockDim.x+threadIdx.x;\n"
            << "for (int i = tid; i < nip; i += blockDim.x*gridDim.x) {\n";
          // s << "for (size_t i = 0; i < nip; i++) {\n";
          
          s << code.body << endl;
          
          // missing: last step nr
          for (int j = 0; j < cf->Dimension(); j++)
            s << "output[i+"<<base_output+j<<"*dist_output] = "
              << Var(compiledcf->Steps().Size()-1, j, cf->Dimensions()).code << ";\n";
          base_output += cf->Dimension();
          
          s << "}\n}";
        }
      s << "}\n";

      s << 
        "extern \"C\" void ApplyIPFunction (size_t nip, double * input, size_t dist_input,\n"
        "                      double * output, size_t dist_output) {\n"
        "  ApplyIPFunctionKernel<<<256,256>>> (nip, input, dist_input, output, dist_output); } \n";

      
      
      cout << IM(9) << s.str() << endl;

      // CUDA - compilation:

      auto dir = CreateTempDir();
      auto prefix = dir.append("GPUcode");
      auto src_file = filesystem::path(prefix).concat(".cu").u8string();
      auto lib_file = filesystem::path(prefix).concat(".so").u8string();
      
      ofstream codefile(src_file);
      codefile << s.str();
      codefile.close();
      cout << "calling compiler, srcfile = " << src_file << ", libfile = " << lib_file << endl;
      return;
      int err = system( ("ngs_nvcc -shared -Xcompiler -fPIC " + src_file + " -o "+lib_file).c_str() );
      if (err) throw Exception ("problem calling compiler");
      system( "ls -alt");
      system( "pwd");
      library = make_unique<SharedLibrary>(lib_file, dir);
      compiled_function = library->GetFunction<lib_function> ("ApplyIPFunction");
    }
    
    virtual void Mult (const BaseVector & x, BaseVector & y) const override
    {
      const UnifiedVector & ux = dynamic_cast<const UnifiedVector&> (x);
      UnifiedVector & uy = dynamic_cast<UnifiedVector&> (y);
      
      ux.UpdateDevice();
      uy.UpdateDevice();

      compiled_function(nip, ux.DevData(), nip, uy.DevData(), nip);
      cudaDeviceSynchronize();

    }

    virtual int VHeight() const override { return h; }
    virtual int VWidth() const { return w; }    
  };


  
  void InitApplyIntegrationPoints ()
  {
    BaseMatrix::RegisterDeviceMatrixCreator(typeid(ApplyIntegrationPoints),
                                            [] (const BaseMatrix & bmat) -> shared_ptr<BaseMatrix>
                                            {
                                              auto & mat = dynamic_cast<const ApplyIntegrationPoints&>(bmat);
                                              return make_shared<DevApplyIntegrationPoints>(mat);
                                            });
  }
};
