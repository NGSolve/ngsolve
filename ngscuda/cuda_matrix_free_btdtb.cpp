
#include <la.hpp>
#include <comp.hpp>
#include <memory>
#include "cuda_linalg.hpp"
#include "cuda_profiler.hpp"

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

    unique_ptr<Matrix<Dev<double>>> dev_bx, dev_by;
    unique_ptr<Array<Dev<int>>> dev_dofx, dev_dofy;
    unique_ptr<Array<Dev<double>>> dev_d, dev_jac;

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
      dev_bx = make_unique<Matrix<Dev<double>>> (compress(mat.Bx, nzxk, nzxc));
      dev_by = make_unique<Matrix<Dev<double>>> (compress(mat.By, nzyk, nzyc));

      Array<int> flatdofx(nel*locdofsx), flatdofy(nel*locdofsy);
      for (auto i : Range(nel))
        {
          for (auto k : Range(locdofsx)) flatdofx[i*locdofsx+k] = mat.dofx[i][k];
          for (auto k : Range(locdofsy)) flatdofy[i*locdofsy+k] = mat.dofy[i][k];
        }
      dev_dofx = make_unique<Array<Dev<int>>> (flatdofx);
      dev_dofy = make_unique<Array<Dev<int>>> (flatdofy);

      dev_d = make_unique<Array<Dev<double>>> (FlatArray<double>(mat.D.GetTotalSize(), mat.D.Data()));
      dev_jac = make_unique<Array<Dev<double>>> (FlatArray<double>(mat.Jacobi.GetTotalSize(), mat.Jacobi.Data()));

      // generate cuda code:

      auto genapply = [](const Array<shared_ptr<DifferentialOperator>> & diffops,
                         string refvar, string worldvar, bool trans)
      {
        stringstream s;
        int starti = 0, startiref = 0;
        for (auto & dop : diffops)
          {
            string refexpr   = "FlatVec<"+ToString(dop->DimRef())+">("+refvar+"+"+ToString(startiref)+")";
            string worldexpr = "FlatVec<"+ToString(dop->Dim())+">("+worldvar+"+"+ToString(starti)+")";
            string invar  = trans ? worldexpr : refexpr;
            string outvar = trans ? refexpr   : worldexpr;
            s << "          " << dop->GenerateTransformationCode(invar, outvar, trans);
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

      // all 32 lanes join every butterfly, lanes with j >= NIP contribute zero;
      // after the xor-reduction every lane holds the sum, lane k%32 writes it
      stringstream sby;
      sby << "          GenRegionTracer rt(brt, 5);\n";
      sby << "          double v_reduced = 0.0;\n";
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
            sby << "          for (int off = 16; off; off >>= 1)\n";
            sby << "            v_reduced += __shfl_xor_sync(0xffffffff, v_reduced, off);\n";
            sby << "          if (threadIdx.x == " << k%32 << ") s_vecy[" << k << "] += v_reduced;\n";
          }
        }

      Code allcode;

      const bool only_loadstore = mat.opts.only_loadstore;
      const bool use_atomic     = mat.opts.atomic;

      stringstream computeblock;
      if (!only_loadstore)
        computeblock <<
          "      Mat<SDIM,SDIM> F;\n"
          "      double J;\n"
          "      {\n"
          "      GenRegionTracer rt(brt, 2);\n"
          "      for (int i = 0; i < SDIM*SDIM; i++)\n"
          "        F(i/SDIM,i%SDIM) = jacs[(size_t(el)*SDIM*SDIM+i)*NIPJ];\n"
          "      J = " << detexpr << ";\n"
          "      }\n"
          "      {\n"
          "      GenRegionTracer rt(brt, 3);\n"
          "      for (int j0 = 0; j0 < NIP; j0 += blockDim.x)\n"
          "        {\n"
          "          int j = j0 + threadIdx.x;\n"
          "          double py[DIMYREF] = {};\n"
          "          if (j < NIP)\n"
          "            {\n"
          "          double px[DIMXREF] = {};\n"
          << sbx.str() <<
          "          double hx[DIMX];\n"
          << genapply(mat.diffopsx, "px", "hx", false) <<
          "          const double * Del = dvals + size_t(el)*DIMY*DIMX*NIPD;\n"
          "          double hy[DIMY];\n"
          "          for (int r = 0; r < DIMY; r++)\n"
          "            {\n"
          "              double sum = 0.0;\n"
          "              for (int c = 0; c < DIMX; c++)\n"
          "                sum += Del[(r*DIMX+c)*NIPD + " << (nipD==1 ? "0" : "j") << "]*hx[c];\n"
          "              hy[r] = weights[j]*sum;\n"
          "            }\n"
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
        "  GenBlockTracer brt;\n"
        "  __shared__ double s_vecx[LOCDOFSX];\n"
        "  __shared__ double s_vecy[LOCDOFSY];\n"
        "  for (int el = blockIdx.x; el < NEL; el += gridDim.x)\n"
        "    {\n"
        "      __syncthreads();\n"
        "      {\n"
        "      GenRegionTracer rt(brt, 1);\n"
        "      for (int i = threadIdx.x; i < LOCDOFSX; i += blockDim.x)\n"
        "        s_vecx[i] = input[dofx[size_t(el)*LOCDOFSX+i]];\n"
        "      for (int i = threadIdx.x; i < LOCDOFSY; i += blockDim.x)\n"
        "        s_vecy[i] = 0.0;\n"
        "      }\n"
        "      __syncthreads();\n"
        << computeblock.str() <<
        "      __syncthreads();\n"
        "      {\n"
        "      GenRegionTracer rt(brt, 4);\n"
        "      for (int i = threadIdx.x; i < LOCDOFSY; i += blockDim.x)\n"
        << scatter <<
        "      }\n"
        "    }\n"
        "}\n"
        "extern \"C\" void MatrixFreeBTDTBFunction (double s, BareVector<Dev<double>> input,\n"
        "                      BareVector<Dev<double>> output, cudaStream_t stream) {\n"
        "  MatrixFreeBTDTBKernel<<<NBLOCKS,32,0,stream>>> (s, (double*)input.Data(), (double*)output.Data()); }\n";

      allcode.body = body.str();

      allcode.AddPointer(dev_bx->Data(), "bxvals", "double *", "__device__");
      allcode.AddPointer(dev_by->Data(), "byvals", "double *", "__device__");
      allcode.AddPointer(dev_dofx->Data(), "dofx", "int *", "__device__");
      allcode.AddPointer(dev_dofy->Data(), "dofy", "int *", "__device__");
      allcode.AddPointer(dev_d->Data(), "dvals", "double *", "__device__");
      allcode.AddPointer(dev_jac->Data(), "jacs", "double *", "__device__");
      #ifdef NGS_CUDA_DEVICE_TIMERS
      allcode.AddPointer(ngs_cuda::GetDevTraceDataPtr(), "tr_trace", "void *", "__device__");
      allcode.AddPointer(ngs_cuda::GetDevTraceBlockDataPtr(), "tr_block", "void *", "__device__");
      allcode.AddPointer(ngs_cuda::GetDevTraceStatePtr(), "tr_state", "void *", "__device__");
      #endif // NGS_CUDA_DEVICE_TIMERS

      stringstream s_top;
      s_top <<
        "#include <bla.hpp>\n"
        "#include <cuda_core.hpp>\n"
        "#include <cuda_linalg.hpp>\n"
        "#include <cuda_profiler.hpp>\n"
        "#include <cstddef>\n"
        "using namespace ngbla;\n";

      // device timers: same tracing as DeviceBlockRegionTracer/DeviceRegionTracer,
      // but writing through raw pointers into the buffers of the main cuda module
      s_top << R"RAW(
#ifdef NGS_CUDA_DEVICE_TIMERS
struct GenBlockTracer
{
  int blockNr, threadNr;
  unsigned trace_pos, trace_end;
  __device__ GenBlockTracer ()
    : blockNr(blockIdx.x), threadNr(threadIdx.x), trace_pos(0), trace_end(0)
  {
    if (threadNr == 0)
      {
        auto & bd = ((ngs_cuda::DevTraceBlockData*)tr_block)[blockNr];
        bd.start = clock64();
        bd.start_global = ngs_cuda::GetGlobalTimer();
        bd.start_smid = ngs_cuda::GetSMID();
        if (blockNr == 0)
          ((ngs_cuda::DevTraceState*)tr_state)->nblocks = gridDim.x;
      }
    __syncwarp();
  }
  __device__ unsigned Alloc ()
  {
    using namespace ngs_cuda;
    if (trace_pos == trace_end)
      {
        trace_pos = TRACE_CHUNK_SIZE * atomicAdd(&((DevTraceState*)tr_state)->chunk_counter, 1u);
        trace_end = trace_pos + TRACE_CHUNK_SIZE;
      }
    if (trace_pos >= N_MAX_TRACER_OBJECTS)
      {
        trace_pos = N_MAX_TRACER_OBJECTS;
        trace_end = N_MAX_TRACER_OBJECTS + 1;
        return N_MAX_TRACER_OBJECTS;
      }
    return trace_pos++;
  }
  __device__ ~GenBlockTracer ()
  {
    if (threadNr == 0)
      {
        auto & bd = ((ngs_cuda::DevTraceBlockData*)tr_block)[blockNr];
        bd.stop = clock64();
        bd.stop_global = ngs_cuda::GetGlobalTimer();
        bd.stop_smid = ngs_cuda::GetSMID();
      }
  }
};
struct GenRegionTracer
{
  bool active;
  unsigned id;
  __device__ GenRegionTracer (GenBlockTracer & tr, int timer_nr)
    : active(tr.threadNr == 0), id(active ? tr.Alloc() : ngs_cuda::N_MAX_TRACER_OBJECTS)
  {
    if (active)
      {
        auto & td = ((ngs_cuda::DevTraceData*)tr_trace)[id];
        td.blockNr = tr.blockNr;
        td.timer_nr = timer_nr;
        td.start = clock64();
      }
  }
  __device__ ~GenRegionTracer ()
  {
    if (active)
      ((ngs_cuda::DevTraceData*)tr_trace)[id].stop = clock64();
  }
};
#else // NGS_CUDA_DEVICE_TIMERS
struct GenBlockTracer { __device__ GenBlockTracer () {} };
struct GenRegionTracer { __device__ GenRegionTracer (GenBlockTracer &, int) {}; };
#endif // NGS_CUDA_DEVICE_TIMERS
)RAW";
      s_top << "static constexpr int NBLOCKS = " << std::min(nel, size_t(4096))
            << ", NEL = " << nel << ", NIP = " << nip << ", NIPD = " << nipD
            << ", LOCDOFSX = " << locdofsx << ", LOCDOFSY = " << locdofsy
            << ", DIMX = " << dimx << ", DIMY = " << dimy
            << ", DIMXREF = " << dimxref << ", DIMYREF = " << dimyref
            << ", SDIM = " << dims << ", NIPJ = " << nipJ << ";\n";
      s_top.precision(17);
      s_top << "__constant__ double weights[NIP] = {";
      for (size_t j = 0; j < nip; j++)
        s_top << (j ? "," : "") << mat.weights(j);
      s_top << "};\n";

      allcode.top += s_top.str();

      cout << IM(9) << allcode.top << allcode.body << endl;

      if (!mat.opts.write_kernel.empty())
        {
          ofstream out(mat.opts.write_kernel);
          out << allcode.top << allcode.body << allcode.pointer;
          cout << IM(3) << "wrote generated GPU kernel to " << mat.opts.write_kernel << endl;
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
      if (synckernels) cudaDeviceSynchronize();
    }

    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override
    {
      static Timer timer("DevMatrixFreeBTDTB::MultAdd");
      UnifiedVectorWrapper ux(x);
      UnifiedVectorWrapper uy(y);

      ngs_cuda::CudaRegionTimer cutimer(timer, &timer_names);
      compiled_function(s, ux.FVDevRO(), uy.FVDev(), ngs_cuda_stream);
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
