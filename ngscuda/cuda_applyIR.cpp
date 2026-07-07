
#include <la.hpp>
#include <comp.hpp>
#include <memory>
#include "cuda_linalg.hpp"

using namespace ngcomp;

namespace ngla
{
  extern bool synckernels;
  
  class DevApplyIntegrationPoints : public DevMatrix
  {
    size_t h, w;
    size_t dimx, dimy, nip;
    unique_ptr<SharedLibrary> library;

    typedef void (*lib_function)(size_t nip, BareVector<Dev<double>> input, size_t dist_input,
                                 BareVector<Dev<double>> output, size_t dist_output, cudaStream_t stream);
    lib_function compiled_function = nullptr;

    unique_ptr<Matrix<Dev<double>>> dev_points, dev_normals, dev_coefs;
    unique_ptr<Array<Dev<int>>> dev_elidx;

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
      auto & input_coefs = aipmat.GetInputCoefs();
      auto & input_coef_offset = aipmat.GetInputCoefOffset();

      Array<int> proxyoffset;
      int starti = 0;
      for (auto proxy : trialproxies)
        {
          proxyoffset.Append (starti);
          starti += proxy->Evaluator()->Dim();
        }

      Code allcode;

      int base_output = 0;
      for (auto cf : aipmat.GetCFs())
        {
          auto compiledcf = Compile (cf, false);
          Code code = compiledcf->GenerateProgram(0, false);
          stringstream s;
          
          s << "{\n";
          // cout << code.header << endl;
          
          for (auto step : Range(compiledcf->Steps()))
            {
              auto stepcf = compiledcf->Steps()[step];

              bool uses_values =
                code.body.find("values_" + ToString(step) + "(") != string::npos;

              if (auto proxycf = dynamic_cast<ProxyFunction*> (stepcf))
                {
                  auto pos = trialproxies.Pos(proxycf);
                  if (pos != trialproxies.ILLEGAL_POSITION)
                    {
                      s << "auto values_" << step << " = [dist_input,input](size_t i, int comp)\n"
                        " { return input[i + (comp+" << proxyoffset[pos] << ")*dist_input]; };\n";
                      s << "bool constexpr has_values_" << step << " = true;\n" << endl;
                      for (int i = 0; i < proxycf->Dimension(); i++)
                        s << Var("comp", step,i,proxycf->Dimensions()).Declare("double", 0.0);
                    }
                }
              else if (uses_values)
                {
                  auto pos = input_coefs.Pos(stepcf);
                  int off = (pos != input_coefs.ILLEGAL_POSITION) ? input_coef_offset[pos] : 0;
                  s << "auto values_" << step << " = [](size_t i, int comp)\n"
                    " { return coef_input[i + (comp+" << off << ")*dist_coef]; };\n" << endl;
                }
            }

          s << "[[maybe_unused]] auto points = [](size_t i, int comp)\n"
            " { return pnts[i+comp*dist_pnts]; };\n";
          s << "[[maybe_unused]] auto normals = [](size_t i, int comp)\n"
            " { return nvs[i+comp*dist_normals]; };\n";
        

          s << "int tid = blockIdx.x*blockDim.x+threadIdx.x;\n"
            << "for (int i = tid; i < nip; i += blockDim.x*gridDim.x) {\n";
          // s << "for (size_t i = 0; i < nip; i++) {\n";

          // domain-wise coefficient functions switch on 'domain_index'
          if (code.body.find("domain_index") != string::npos)
            s << "[[maybe_unused]] int domain_index = elidx[i];\n";

          s << code.body << endl;
          
          // missing: last step nr
          for (int j = 0; j < cf->Dimension(); j++)
            s << "output[i+"<<base_output+j<<"*dist_output] = "
              << Var(compiledcf->Steps().Size()-1, j, cf->Dimensions()).code << ";\n";
          base_output += cf->Dimension();
          
          s << "}\n}";

          allcode.body += s.str();
        }

      allcode.body += "\n}\n"; // end kernel
      allcode.body +=
        "extern \"C\" void ApplyIPFunction (size_t nip, BareVector<Dev<double>> input, size_t dist_input,\n"
        "                      BareVector<Dev<double>> output, size_t dist_output, cudaStream_t stream) {\n"
        "  ApplyIPFunctionKernel<<<256,256,0,stream>>> (nip, (double*)input.Data(), dist_input, (double*)output.Data(), dist_output); } \n";

      stringstream s_top;
      s_top <<
        "#include <bla.hpp>\n"
        "#include <cuda_core.hpp>\n"        
        "#include <cuda_linalg.hpp>\n"
        "#include <cstddef>\n"
        "using namespace ngbla;\n"
      ;

      auto points = aipmat.GetPoints();
      if(allcode.body.find("points") != string::npos)
        dev_points = make_unique<Matrix<Dev<double>>>(points);
      allcode.AddPointer(dev_points ?  dev_points->Data() : nullptr, "pnts", "double *", "__device__");
      allcode.AddPointer((void*)points.Dist(), "dist_pnts", "uint", "__device__");

      auto normals = aipmat.GetNormals();
      // s_top << "size_t dist_normals = " << normals.Dist() << ";\n";
      if(allcode.body.find("normals") != string::npos)
        dev_normals = make_unique<Matrix<Dev<double>>>(normals);
      allcode.AddPointer(dev_normals ?  dev_normals->Data() : nullptr, "nvs", "double *", "__device__");
      allcode.AddPointer((void*)normals.Dist(), "dist_normals", "uint", "__device__");

      auto coef_values = aipmat.GetCoefValues();
      if(allcode.body.find("coef_input") != string::npos)
        dev_coefs = make_unique<Matrix<Dev<double>>>(coef_values);
      allcode.AddPointer(dev_coefs ?  dev_coefs->Data() : nullptr, "coef_input", "double *", "__device__");
      allcode.AddPointer((void*)coef_values.Dist(), "dist_coef", "uint", "__device__");

      auto element_index = aipmat.GetElementIndex();
      if(allcode.body.find("elidx") != string::npos)
        dev_elidx = make_unique<Array<Dev<int>>>(element_index);
      allcode.AddPointer(dev_elidx ?  dev_elidx->Data() : nullptr, "elidx", "int *", "__device__");

      s_top << "__global__ void ApplyIPFunctionKernel (size_t nip, double * input, size_t dist_input,\n"
        "                      double * output, size_t dist_output) {\n";
      
      allcode.top += s_top.str();
      
      cout << IM(9) << allcode.body << endl;

      // CUDA - compilation:

      auto dir = CreateTempDir();
      auto prefix = dir.append("GPUcode");

      auto src_file = filesystem::path(prefix).concat(".cu");
      auto ptr_file = filesystem::path(prefix).concat("_ptrs.cu");
      
      ofstream{src_file} << allcode.top << allcode.body;
      ofstream{ptr_file} << allcode.pointer;

      library = CompileCode( {src_file, ptr_file}, {}, false, "ngs_nvcc", "ngs_nvlink" );
      compiled_function = library->GetSymbol<lib_function> ("ApplyIPFunction");
    }
    
    virtual void Mult (const BaseVector & x, BaseVector & y) const override
    {
      UnifiedVectorWrapper ux(x);
      UnifiedVectorWrapper uy(y);
      
      // const UnifiedVector & ux = dynamic_cast<const UnifiedVector&> (x);
      // UnifiedVector & uy = dynamic_cast<UnifiedVector&> (y);
      
      // ux.UpdateDevice();
      // uy.UpdateDevice();

      compiled_function(nip, ux.FVDevRO(), nip, uy.FVDev(), nip, ngs_cuda_stream);
      if (synckernels) cudaDeviceSynchronize();
    }

    virtual int VHeight() const override { return h; }
    virtual int VWidth() const override { return w; }
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
