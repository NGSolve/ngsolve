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
    gc->Mult(*r, *z);   // z = M^{-1} r
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

    // ── CUDA graph capture of CG iteration body (excl. convergence check) ──
    CudaGraph g_cg;
    bool use_graph = false;
    if (!getenv("NO_CUDA_GRAPH")) {
        g_cg.BeginCapture();
            ga->Mult(*p, *q);                          // q = A*p
            (*p).InnerProduct(*q, *pq);                // pq = <p,q>
            ualpha     = Div(Scal(urz), Scal(upq));    // alpha = rz/pq
            uneg_alpha = Neg(Scal(ualpha));             // neg_alpha = -alpha
            (*x).Add(*alpha,     *p);                  // x += alpha*p
            (*r).Add(*neg_alpha, *q);                  // r -= alpha*q
            gc->Mult(*r, *z);                          // z = M^{-1}*r
            (*r).InnerProduct(*z, *rz_new);            // rz_new = <r,z>
            ubeta = Div(Scal(urz_new), Scal(urz));     // beta = rz_new/rz
            urz   = Scal(urz_new);                     // rz = rz_new
            (*p).Scale(*beta);                         // p *= beta
            (*p).Add(1.0, *z);                         // p += z
        g_cg.EndCapture();
        SyncNGSStream();
        use_graph = g_cg.IsValid();
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

    sol = *x;
}

} // namespace ngla
