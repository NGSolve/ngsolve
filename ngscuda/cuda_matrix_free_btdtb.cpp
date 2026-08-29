
#include <la.hpp>
#include <comp.hpp>
#include <memory>
#include "cuda_linalg.hpp"
#include "cuda_profiler.hpp"

#include "tinybla.hpp"

using namespace ngcomp;

namespace ngla
{
  extern bool synckernels;

  class DevMatrixFreeBTDTB : public DevMatrix
  {
    size_t h, w;
    shared_ptr<SharedLibrary> library;

    typedef void (*lib_function)(double s, BareVector<Dev<double>> input, BareVector<Dev<double>> output,
                                 cudaStream_t stream);
    lib_function compiled_function = nullptr;

    unique_ptr<Array<Dev<int>>> dev_dofx, dev_dofy;
    Array<unique_ptr<Array<Dev<double>>>> dev_data;
    Array<unique_ptr<Array<Dev<float>>>> dev_data_f;

    // names for the device timer regions in the generated kernel
    static inline std::initializer_list<const char*> timer_names =
      { "load B", "gather x", "trafo", "IP loop", "scatter y", "atomic_add" };

  public:
    DevMatrixFreeBTDTB (const MatrixFreeBTDTB & mat)
    {
      auto [locdofsx, dimxref, nip] = mat.Bx.Shape();
      auto [locdofsy, dimyref, nipy] = mat.By.Shape();
      auto [numels, dimy, dimx, nipD] = mat.D.Shape();
      auto [numelsJ, dimr, dims, nipJ] = mat.Jacobi.Shape();
      size_t nel = mat.elnums.Size();

      if (dimr != dims)
        throw Exception("DevMatrixFreeBTDTB: only volume (square) jacobians supported");

      h = mat.height;
      w = mat.width;

      // compress Bx/By: sparsity pattern is the same for all integration points,
      // store only non-zeros, packed ip-wise: vals(nz, ip)
      Array<int> nzxk, nzxc, nzyk, nzyc;
      auto compress = [](const Tensor<3> & B, Array<int> & nzk, Array<int> & nzc) -> Matrix<>
      {
        auto [nd, dim, np] = B.Shape();
        for (size_t k = 0; k < nd; k++)
          for (size_t c = 0; c < dim; c++)
            for (size_t j = 0; j < np; j++)
              if (B(k,c,j) != 0.0)
                {
                  nzk.Append(k);
                  nzc.Append(c);
                  break;
                }
        Matrix<> vals(nzk.Size(), np);
        for (auto i : Range(nzk))
          for (size_t j = 0; j < np; j++)
            vals(i,j) = B(nzk[i], nzc[i], j);
        return vals;
      };
      const bool fp32      = mat.opts.fp32;
      const bool nonlinear = mat.opts.nonlinear;

      auto upload = [&](FlatArray<double> vals) -> const void*
      {
        if (fp32)
          {
            Array<float> fvals(vals.Size());
            for (auto i : Range(vals)) fvals[i] = float(vals[i]);
            dev_data_f.Append (make_unique<Array<Dev<float>>> (fvals));
            return dev_data_f.Last()->Data();
          }
        dev_data.Append (make_unique<Array<Dev<double>>> (vals));
        return dev_data.Last()->Data();
      };

      Matrix<> bxvals = compress(mat.Bx, nzxk, nzxc);
      const void * bxptr = upload(FlatArray<double>(bxvals.Height()*bxvals.Width(), bxvals.Data()));
      // fold quadrature weights into By, saves the per-lane weights[j] lookup in the kernel
      Matrix<> byvals = compress(mat.By, nzyk, nzyc);
      for (auto i : Range(byvals.Height()))
        for (size_t j = 0; j < nip; j++)
          byvals(i,j) *= mat.weights(j);
      const void * byptr = upload(FlatArray<double>(byvals.Height()*byvals.Width(), byvals.Data()));

      Array<int> flatdofx(nel*locdofsx), flatdofy(nel*locdofsy);
      for (auto i : Range(nel))
        {
          for (auto k : Range(locdofsx)) flatdofx[i*locdofsx+k] = mat.dofx[i][k];
          for (auto k : Range(locdofsy)) flatdofy[i*locdofsy+k] = mat.dofy[i][k];
        }
      dev_dofx = make_unique<Array<Dev<int>>> (flatdofx);
      dev_dofy = make_unique<Array<Dev<int>>> (flatdofy);

      const void * dptr = nullptr;
      if (!nonlinear)
        dptr = upload(FlatArray<double>(mat.D.GetTotalSize(), mat.D.Data()));
      const void * jacptr = upload(FlatArray<double>(mat.Jacobi.GetTotalSize(), mat.Jacobi.Data()));

      // generate cuda code:

      auto genapply = [](const Array<shared_ptr<DifferentialOperator>> & diffops,
                         string refvar, string worldvar, bool trans)
      {
        stringstream s;
        int starti = 0, startiref = 0;
        for (auto & dop : diffops)
          {
            string refexpr   = refvar+"+"+ToString(startiref);
            string worldexpr = worldvar+"+"+ToString(starti);
            string inptr  = trans ? worldexpr : refexpr;
            string outptr = trans ? refexpr   : worldexpr;
            int indim  = trans ? dop->Dim()    : dop->DimRef();
            int outdim = trans ? dop->DimRef() : dop->Dim();
            s << "          {\n";
            s << "          tinybla::Vec<" << outdim << ",SCAL> res;\n";
            s << "          " << dop->GenerateTransformationCode("LoadVec<"+ToString(indim)+">("+inptr+")", "res", trans);
            s << "          StoreVec(" << outptr << ", res);\n";
            s << "          }\n";
            starti    += dop->Dim();
            startiref += dop->DimRef();
          }
        return s.str();
      };

      string detexpr;
      if (dims == 2)
        detexpr = "F(0,0)*F(1,1) - F(0,1)*F(1,0)";
      else
        detexpr = "F(0,0)*(F(1,1)*F(2,2)-F(1,2)*F(2,1)) - F(0,1)*(F(1,0)*F(2,2)-F(1,2)*F(2,0))"
                  " + F(0,2)*(F(1,0)*F(2,1)-F(1,1)*F(2,0))";

      // explicitly unrolled apply of Bx^T / By, indices hard-coded
      stringstream sbx;
      for (auto i : Range(nzxk))
        sbx << "          px[" << nzxc[i] << "] += bxvals[" << i*nip << "+j]*s_vecx[" << nzxk[i] << "];\n";

      stringstream sby;
      sby << "          DeviceRegionTracer rt(brt, 5);\n";
      for (size_t k = 0; k < locdofsy; k++)
        {
          stringstream terms;
          bool first = true;
          for (auto i : Range(nzyk))
            if (nzyk[i] == int(k))
              {
                terms << (first ? "" : " + ") << "byvals[" << i*nip << "+j]*py[" << nzyc[i] << "]";
                first = false;
              }
          if (!first)
          {
            sby << "          v_reduced = (j < NIP) ? " << terms.str() << " : 0.0;\n";
            sby << "          for (int off = BS_IPTS/2; off; off >>= 1)\n";
            sby << "            v_reduced += __shfl_xor_sync(submask, v_reduced, off);\n";
            sby << "          if (j < NIP) s_vecy[" << k << "] += v_reduced;\n";
          }
        }

      string physcode;
      if (nonlinear)
      {
        stringstream phys;
        phys << "          { // physics\n"
                "          [[maybe_unused]] int i = 0;\n";
        for (auto k : Range(mat.test_proxies))
          {
            CoefficientFunction::T_DJC cache;
            auto diffcf = mat.cf -> DiffJacobi (mat.test_proxies[k], cache);
            auto compiledcf = Compile (diffcf, false);
            Code code = compiledcf->GenerateProgram(0, false);

            phys << "          { // equation " << k << "\n";
            for (auto step : Range(compiledcf->Steps()))
              {
                auto stepcf = compiledcf->Steps()[step];
                if (auto proxycf = dynamic_cast<ProxyFunction*> (stepcf))
                  {
                    auto pos = mat.trial_proxies.Pos(proxycf);
                    if (pos != mat.trial_proxies.ILLEGAL_POSITION)
                      {
                        phys << "auto values_" << step << " = [&](int ip, int nr) { return hx[" << mat.ranges_x[pos].First() << "+nr]; };\n";
                        phys << "constexpr bool has_values_" << step << " = true;\n";
                        for (auto l : Range(proxycf->Dimension()))
                          phys << Var("comp", step, l, proxycf->Dimensions()).Declare("double", 0.0);
                      }
                  }
              }
            phys << code.body;

            int steps = compiledcf->Steps().Size();
            IntRange rangey = mat.ranges_y[k];
            if (rangey.Size()==1)
              phys << "hy[" << rangey[0] << "] = var_" << steps-1 << ";\n";
            else
              for (auto l : Range(rangey))
                phys << "hy[" << rangey[l] << "] = var_" << steps-1 << "_" << l << ";\n";
            phys << "          }\n";
          }
        phys << "          }\n";
        physcode = phys.str();
        if (fp32)
          for (size_t p = physcode.find("double"); p != string::npos; p = physcode.find("double", p+5))
            physcode.replace(p, 6, "float");
      }

      Code allcode;

      const bool only_loadstore = mat.opts.only_loadstore;
      const bool use_atomic     = mat.opts.atomic;
      size_t epb                = std::max(1, mat.opts.BS_els);
      
      if (epb * mat.opts.BS_ipts != 32)
          epb = 32 / mat.opts.BS_ipts;

      stringstream dapply;
      if (nonlinear)
        dapply << physcode
               << "          for (int r = 0; r < DIMY; r++) hy[r] *= fabs(J);\n";
      else
        dapply <<
          "          const SCAL * Del = dvals + size_t(el)*DIMY*DIMX*NIPD;\n"
          "          for (int r = 0; r < DIMY; r++)\n"
          "            {\n"
          "              SCAL sum = 0.0;\n"
          "              for (int c = 0; c < DIMX; c++)\n"
          "                sum += Del[(r*DIMX+c)*NIPD + " << (nipD==1 ? "0" : "j") << "]*hx[c];\n"
          "              hy[r] = sum;\n"
          "            }\n";

      stringstream computeblock;
      if (!only_loadstore)
        computeblock <<
          "      tinybla::Mat<SDIM,SDIM,SCAL> F;\n"
          "      [[maybe_unused]] SCAL J;\n"
          "      {\n"
          "      DeviceRegionTracer rt(brt, 2);\n"
          "      for (int i = 0; i < SDIM*SDIM; i++)\n"
          "        F(i/SDIM,i%SDIM) = jacs[(size_t(el)*SDIM*SDIM+i)*NIPJ];\n"
          "      J = " << detexpr << ";\n"
          "      }\n"
          "      {\n"
          "      DeviceRegionTracer rt(brt, 3);\n"
          "      SCAL v_reduced;\n"
          "      for (int j0 = 0; j0 < NIP; j0 += blockDim.x)\n"
          "        {\n"
          "          int j = j0 + threadIdx.x;\n"
          "          SCAL py[DIMYREF] = {};\n"
          "          if (j < NIP)\n"
          "            {\n"
          "          SCAL px[DIMXREF] = {};\n"
          << sbx.str() <<
          "          SCAL hx[DIMX];\n"
          << genapply(mat.diffopsx, "px", "hx", false) <<
          "          SCAL hy[DIMY];\n"
          << dapply.str()
          << genapply(mat.diffopsy, "py", "hy", true) <<
          "            }\n"
          << sby.str() <<
          "        }\n"
          "      }\n";

      const char * scatter = use_atomic
        ? "        atomicAdd(&output[dofy[size_t(el)*LOCDOFSY+i]], s*s_vecy[i]);\n"
        : "        output[dofy[size_t(el)*LOCDOFSY+i]] += s*s_vecy[i];\n";

      stringstream body;
      body <<
        "__global__ void MatrixFreeBTDTBKernel (double s, double * input, double * output)\n"
        "{\n"
        "  DeviceBlockRegionTracer brt(gridDim, blockDim, blockIdx, threadIdx);\n"
        "  __shared__ SCAL s_vecx_all[BS_ELS*LOCDOFSX];\n"
        "  __shared__ SCAL s_vecy_all[BS_ELS*LOCDOFSY];\n"
        "  const int ey = threadIdx.y;\n"
        "  SCAL * s_vecx = s_vecx_all + ey*LOCDOFSX;\n"
        "  SCAL * s_vecy = s_vecy_all + ey*LOCDOFSY;\n"
        "  // sync only the BS_IPTS threads working on this element\n"
        "  const unsigned submask = unsigned((1ull << BS_IPTS) - 1) << (BS_IPTS*ey);\n"
        "  for (int el = blockIdx.x*BS_ELS + ey; el < NEL; el += gridDim.x*BS_ELS)\n"
        "    {\n"
        "      __syncwarp(submask);\n"
        "      {\n"
        "      DeviceRegionTracer rt(brt, 1);\n"
        "      for (int i = threadIdx.x; i < LOCDOFSX; i += blockDim.x)\n"
        "        s_vecx[i] = SCAL(input[dofx[size_t(el)*LOCDOFSX+i]]);\n"
        "      for (int i = threadIdx.x; i < LOCDOFSY; i += blockDim.x)\n"
        "        s_vecy[i] = 0.0;\n"
        "      }\n"
        "      __syncwarp(submask);\n"
        << computeblock.str() <<
        "      __syncwarp(submask);\n"
        "      {\n"
        "      DeviceRegionTracer rt(brt, 4);\n"
        "      for (int i = threadIdx.x; i < LOCDOFSY; i += blockDim.x)\n"
        << scatter <<
        "      }\n"
        "    }\n"
        "}\n"
        "extern \"C\" void MatrixFreeBTDTBFunction (double s, BareVector<Dev<double>> input,\n"
        "                      BareVector<Dev<double>> output, cudaStream_t stream) {\n"
        "  MatrixFreeBTDTBKernel<<<NBLOCKS,dim3(BS_IPTS,BS_ELS),0,stream>>> (s, (double*)input.Data(), (double*)output.Data()); }\n";

      allcode.body = body.str();


      stringstream s_top;

      if(CudaRegionTimer::IsTimingEnabled())
          s_top << "#define NGS_CUDA_DEVICE_TIMERS\n";

      s_top <<
        "#include <bla.hpp>\n"
        "#include <cuda_core.hpp>\n"
        "#include <cuda_linalg.hpp>\n"
        "#include <cuda_profiler.hpp>\n"
        "#include <cstddef>\n"
        "using namespace ngbla;\n"
        "using namespace ngs_cuda;\n";

      s_top << code_tinybla <<
        "using tinybla::ToMat;\n"
        "template <int S, typename T> __host__ __device__ tinybla::Vec<S,T> LoadVec (const T * p)\n"
        "{ tinybla::Vec<S,T> r; for (int i = 0; i < S; i++) r(i) = p[i]; return r; }\n"
        "template <int S, typename T> __host__ __device__ void StoreVec (T * p, tinybla::Vec<S,T> v)\n"
        "{ for (int i = 0; i < S; i++) p[i] = v(i); }\n"
        "typedef " << (fp32 ? "float" : "double") << " SCAL;\n";

      const size_t nblocks = std::min((nel+epb-1)/epb, size_t(4096));
      s_top << "static constexpr int NBLOCKS = " << nblocks
            << ", NEL = " << nel << ", NIP = " << nip
            << (nonlinear ? string() : ", NIPD = "+ToString(nipD))
            << ", LOCDOFSX = " << locdofsx << ", LOCDOFSY = " << locdofsy
            << ", DIMX = " << dimx << ", DIMY = " << dimy
            << ", DIMXREF = " << dimxref << ", DIMYREF = " << dimyref
            << ", SDIM = " << dims << ", NIPJ = " << nipJ
            << ", BS_ELS = " << epb << ", BS_IPTS = " << mat.opts.BS_ipts
            << ";\n";
      allcode.top += s_top.str();

      const string scalptr = fp32 ? "float *" : "double *";
      allcode.AddPointer(bxptr, "bxvals", scalptr, "__device__");
      allcode.AddPointer(byptr, "byvals", scalptr, "__device__");
      allcode.AddPointer(dev_dofx->Data(), "dofx", "int *", "__device__");
      allcode.AddPointer(dev_dofy->Data(), "dofy", "int *", "__device__");
      if (dptr)
        allcode.AddPointer(dptr, "dvals", scalptr, "__device__");
      allcode.AddPointer(jacptr, "jacs", scalptr, "__device__");
      allcode.pointer += "struct DevTimerData;\nstruct DevTraceData;\nstruct DevTraceBlockData;\nstruct DevTraceState;\n";
      allcode.AddPointer(ngs_cuda::GetDevTimerDataPtr(), "d_timer_data", "DevTimerData *", "__device__");
      allcode.AddPointer(ngs_cuda::GetDevTraceDataPtr(), "d_trace_data", "DevTraceData *", "__device__");
      allcode.AddPointer(ngs_cuda::GetDevTraceBlockDataPtr(), "d_block_data", "DevTraceBlockData *", "__device__");
      allcode.AddPointer(ngs_cuda::GetDevTraceStatePtr(), "d_trace_state", "DevTraceState *", "__device__");


      cout << IM(9) << allcode.top << allcode.body << endl;

      if (mat.opts.write_GPU_kernel)
        {
          ofstream out(*mat.opts.write_GPU_kernel);
          out << allcode.top << allcode.body << allcode.pointer;
          cout << IM(3) << "wrote generated GPU kernel to " << *mat.opts.write_GPU_kernel << endl;
        }

      // CUDA - compilation:

      auto dir = CreateTempDir();
      auto prefix = dir.append("GPUcode");

      auto src_file = filesystem::path(prefix).concat(".cu");
      auto ptr_file = filesystem::path(prefix).concat("_ptrs.cu");

      ofstream{src_file} << allcode.top << allcode.body;
      ofstream{ptr_file} << allcode.pointer;

      library = CompileCode( {src_file, ptr_file}, {}, false, "ngs_nvcc", "ngs_nvlink" );
      compiled_function = library->GetSymbol<lib_function> ("MatrixFreeBTDTBFunction");
    }

    virtual void Mult (const BaseVector & x, BaseVector & y) const override
    {
      static Timer timer("DevMatrixFreeBTDTB::Mult");
      UnifiedVectorWrapper ux(x);
      UnifiedVectorWrapper uy(y);

      uy = 0.0;
      ngs_cuda::CudaRegionTimer cutimer(timer, &timer_names);
      compiled_function(1.0, ux.FVDevRO(), uy.FVDev(), ngs_cuda_stream);
      if (auto err = cudaGetLastError(); err != cudaSuccess)
        throw Exception(string("MatrixFreeBTDTB kernel launch failed: ")+cudaGetErrorString(err));
      if (synckernels) cudaDeviceSynchronize();
    }

    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override
    {
      static Timer timer("DevMatrixFreeBTDTB::MultAdd");
      UnifiedVectorWrapper ux(x);
      UnifiedVectorWrapper uy(y);

      ngs_cuda::CudaRegionTimer cutimer(timer, &timer_names);
      compiled_function(s, ux.FVDevRO(), uy.FVDev(), ngs_cuda_stream);
      if (auto err = cudaGetLastError(); err != cudaSuccess)
        throw Exception(string("MatrixFreeBTDTB kernel launch failed: ")+cudaGetErrorString(err));
      if (synckernels) cudaDeviceSynchronize();
    }

    virtual int VHeight() const override { return h; }
    virtual int VWidth() const override { return w; }
  };



  void InitBTDTB ()
  {
    BaseMatrix::RegisterDeviceMatrixCreator(typeid(MatrixFreeBTDTB),
                                            [] (const BaseMatrix & bmat) -> shared_ptr<BaseMatrix>
                                            {
                                              auto & mat = dynamic_cast<const MatrixFreeBTDTB&>(bmat);
                                              return make_shared<DevMatrixFreeBTDTB>(mat);
                                            });
  }
};
