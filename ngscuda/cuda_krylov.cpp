#include "cuda_linalg.hpp"
#include "cuda_core.hpp"
#include <core/ngcore.hpp>
#include <cmath>

using namespace ngcore;

extern cudaStream_t ngs_cuda_stream;

namespace ngla
{

static void SyncNGSStream()
{
    cudaStreamSynchronize(ngs_cuda::ngs_cuda_stream);
}

void DevCGSolver::Mult(const BaseVector& rhs, BaseVector& sol) const
{
    auto ga = a_dev;
    auto gc = c_dev;
    const int ndof = rhs.Size();

    // pre-allocate working vectors as UnifiedVector
    auto p = make_unique<UnifiedVector>(ndof);
    auto q = make_unique<UnifiedVector>(ndof);
    auto r = make_unique<UnifiedVector>(ndof);
    auto z = make_unique<UnifiedVector>(ndof);
    auto x = make_unique<UnifiedVector>(ndof);

    // pre-allocate GPU-resident scalars
    auto rz        = (*r).CreateScalar();
    auto pq        = (*r).CreateScalar();
    auto rz_new    = (*r).CreateScalar();
    auto alpha     = (*r).CreateScalar();
    auto neg_alpha = (*r).CreateScalar();
    auto beta      = (*r).CreateScalar();

    // cast to UnifiedScalar for expression templates
    auto& urz        = dynamic_cast<UnifiedScalar&>(*rz);
    auto& upq        = dynamic_cast<UnifiedScalar&>(*pq);
    auto& urz_new    = dynamic_cast<UnifiedScalar&>(*rz_new);
    auto& ualpha     = dynamic_cast<UnifiedScalar&>(*alpha);
    auto& uneg_alpha = dynamic_cast<UnifiedScalar&>(*neg_alpha);
    auto& ubeta      = dynamic_cast<UnifiedScalar&>(*beta);

    // ── initialise ──────────────────────────────────────────
    *x = 0.0;
    *r = rhs;
    gc->Mult(*r, *z);
    *p = *z;
    (*r).InnerProduct(*z, *rz);
    SyncNGSStream();

    double r0norm = sqrt(std::abs(rz->GetD()));
    if (r0norm == 0.0) { sol = *x; return; }

    // warm-up: force lazy CUDA init before capture
    ga->Mult(*p, *q);
    SyncNGSStream();

    // ensure all vectors are on device before capture
    dynamic_cast<UnifiedVector&>(*p).UpdateDevice();
    dynamic_cast<UnifiedVector&>(*q).UpdateDevice();
    dynamic_cast<UnifiedVector&>(*r).UpdateDevice();
    dynamic_cast<UnifiedVector&>(*z).UpdateDevice();
    dynamic_cast<UnifiedVector&>(*x).UpdateDevice();
    SyncNGSStream();

    bool use_while_graph = false;
    bool use_graph = false;
    ngs_cuda::CudaWhileGraph g_while;
    CudaGraph g_cg_for_while;
    CudaGraph g_cg;
    int* iter_count_dev = nullptr;

    if (!getenv("NO_CUDA_GRAPH")) {
        double tol = GetPrecision() * r0norm;

        // ── try cudaGraphCondTypeWhile capture ──
        try {
            // First capture iteration body into regular graph (cuSPARSE works here)
            g_cg_for_while.BeginCapture();
                cublasSetPointerMode(Get_CuBlas_Handle(), CUBLAS_POINTER_MODE_DEVICE);
                ga->Mult(*p, *q);
                (*p).InnerProduct(*q, *pq);
                ualpha     = Div(Scal(urz), Scal(upq));
                uneg_alpha = Neg(Scal(ualpha));
                (*x).Add(*alpha,     *p);
                (*r).Add(*neg_alpha, *q);
                gc->Mult(*r, *z);
                (*r).InnerProduct(*z, *rz_new);
                ubeta = Div(Scal(urz_new), Scal(urz));
                urz   = Scal(urz_new);
                (*p).Scale(*beta);
                (*p).Add(1.0, *z);
            g_cg_for_while.EndCapture();
            SyncNGSStream();
            cublasSetPointerMode(Get_CuBlas_Handle(), CUBLAS_POINTER_MODE_HOST);
            // Allocate GPU-side iteration counter
            if (!iter_count_dev) cudaMalloc(&iter_count_dev, sizeof(int));
            cudaMemset(iter_count_dev, 0, sizeof(int));
            // Then build WHILE graph using captured graph as child node
            g_while.Build(g_cg_for_while.GetGraph(), urz.DevPtr(), tol, iter_count_dev, maxsteps);
            SyncNGSStream();
            use_while_graph = g_while.IsValid();
        } catch (ngstd::Exception& e) {
            std::cerr << "[DevCGSolver] CudaWhileGraph capture failed: "
                      << e.What() << std::endl;
            std::cerr << "[DevCGSolver] falling back to iteration graph" << std::endl;
        } catch (...) {
            std::cerr << "[DevCGSolver] CudaWhileGraph capture failed (unknown exception)" << std::endl;
        }

        // ── fall back to per-iteration graph ──
        if (!use_while_graph) {
            g_cg.BeginCapture();
                ga->Mult(*p, *q);
                (*p).InnerProduct(*q, *pq);
                ualpha     = Div(Scal(urz), Scal(upq));
                uneg_alpha = Neg(Scal(ualpha));
                (*x).Add(*alpha,     *p);
                (*r).Add(*neg_alpha, *q);
                gc->Mult(*r, *z);
                (*r).InnerProduct(*z, *rz_new);
                ubeta = Div(Scal(urz_new), Scal(urz));
                urz   = Scal(urz_new);
                (*p).Scale(*beta);
                (*p).Add(1.0, *z);
            g_cg.EndCapture();
            SyncNGSStream();
            use_graph = g_cg.IsValid();
        }
    }

    // re-initialise after warm-up
    *x = 0.0;
    *r = rhs;
    gc->Mult(*r, *z);
    *p = *z;
    (*r).InnerProduct(*z, *rz);
    SyncNGSStream();

    // ── main loop ────────────────────────────────────────────
    int step = 0;

    if (use_while_graph) {
        // entire CG solve on GPU — single graph launch
        cudaMemset(iter_count_dev, 0, sizeof(int));
        g_while.Launch();
        SyncNGSStream();
        cublasSetPointerMode(Get_CuBlas_Handle(), CUBLAS_POINTER_MODE_HOST);
        // read back iteration count for GetSteps()
        cudaMemcpy(&step, iter_count_dev, sizeof(int), cudaMemcpyDeviceToHost);
        cudaFree(iter_count_dev);
        steps = step;
        sol = *x;
        return;
    } else {
        for (int iter = 0; iter < maxsteps; iter++)
        {
            if (use_graph) {
                cublasSetPointerMode(Get_CuBlas_Handle(), CUBLAS_POINTER_MODE_DEVICE);
                g_cg.Launch();
                SyncNGSStream();
                cublasSetPointerMode(Get_CuBlas_Handle(), CUBLAS_POINTER_MODE_HOST);
            } else {
                ga->Mult(*p, *q);
                (*p).InnerProduct(*q, *pq);
                ualpha     = Div(Scal(urz),    Scal(upq));
                uneg_alpha = Neg(Scal(ualpha));
                (*x).Add(*alpha,     *p);
                (*r).Add(*neg_alpha, *q);
                gc->Mult(*r, *z);
                (*r).InnerProduct(*z, *rz_new);
                ubeta = Div(Scal(urz_new), Scal(urz));
                urz   = Scal(urz_new);
                (*p).Scale(*beta);
                (*p).Add(1.0, *z);
                SyncNGSStream();
            }

            step++;
            double res = sqrt(std::abs(rz->GetD()));

            if (printrates)
                cout << "CG iter " << step << "  res = " << res << endl;

            if (res <= GetPrecision() * r0norm || step >= GetMaxSteps()) {
                sol = *x;
                if (printrates)
                    cout << "CG " << (res <= GetPrecision()*r0norm ? "converged" : "max iters")
                         << " after " << step << " iters"
                         << "  res/r0 = " << res/r0norm << endl;
                return;
            }
        }
    }

    steps = step;
    sol = *x;
}

} // namespace ngla
