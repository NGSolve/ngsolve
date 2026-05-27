#include "cuda_linalg.hpp"
#include "cuda_core.hpp"
#include "linalg_kernels.hpp"
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

    // initialise
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
        // diamond + WHILE graph 
        if (getenv("USE_DIAMOND_GRAPH")) {
            try {
                CudaGraph g_pre, g_x, g_r, g_post;

                g_pre.BeginCapture();
                    cublasSetPointerMode(Get_CuBlas_Handle(), CUBLAS_POINTER_MODE_DEVICE);
                    ga->Mult(*p, *q);
                    (*p).InnerProduct(*q, *pq);
                    ualpha     = Div(Scal(urz), Scal(upq));
                    uneg_alpha = Neg(Scal(ualpha));
                g_pre.EndCaptureOnly();

                g_x.BeginCapture();
                    (*x).Add(*alpha, *p);
                g_x.EndCaptureOnly();

                g_r.BeginCapture();
                    (*r).Add(*neg_alpha, *q);
                g_r.EndCaptureOnly();

                g_post.BeginCapture();
                    gc->Mult(*r, *z);
                    (*r).InnerProduct(*z, *rz_new);
                    ubeta = Div(Scal(urz_new), Scal(urz));
                    urz   = Scal(urz_new);
                    (*p).Scale(*beta);
                    (*p).Add(1.0, *z);
                g_post.EndCaptureOnly();

                SyncNGSStream();
                cublasSetPointerMode(Get_CuBlas_Handle(), CUBLAS_POINTER_MODE_HOST);

                // Build diamond from 4 sub-graphs
                ngs_cuda::CudaDiamondGraph g_diamond;
                cudaGraph_t branches[] = {g_x.GetGraph(), g_r.GetGraph()};
                g_diamond.Build(g_pre.GetGraph(),
                                ngcore::FlatArray<cudaGraph_t>(2, branches),
                                g_post.GetGraph());

                // Pass diamond as iteration body to WHILE graph
                if (!iter_count_dev) cudaMalloc(&iter_count_dev, sizeof(int));
                cudaMemset(iter_count_dev, 0, sizeof(int));
                g_while.Build(g_diamond.GetGraph(), urz.DevPtr(), tol, iter_count_dev, maxsteps);
                SyncNGSStream();
                use_while_graph = g_while.IsValid();
            } catch (ngstd::Exception& e) {
                std::cerr << "[DevCGSolver] Diamond capture failed: "
                          << e.What() << std::endl;
                std::cerr << "[DevCGSolver] falling back to WHILE graph" << std::endl;
            } catch (...) {
                std::cerr << "[DevCGSolver] Diamond capture failed (unknown exception)" << std::endl;
            }
        }


        if (!use_while_graph)
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

        // fall back to per-iteration graph
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

    // main loop 
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

void DevTFQMRSolver::Mult(const BaseVector& rhs, BaseVector& sol) const
{
    auto ga = a_dev;
    auto gc = c_dev;
    const int ndof = rhs.Size();

    auto r    = make_unique<UnifiedVector>(ndof);
    auto u    = make_unique<UnifiedVector>(ndof);
    auto v    = make_unique<UnifiedVector>(ndof);
    auto w    = make_unique<UnifiedVector>(ndof);
    auto uhat = make_unique<UnifiedVector>(ndof);
    auto d    = make_unique<UnifiedVector>(ndof);
    auto x    = make_unique<UnifiedVector>(ndof);
    auto tmp  = make_unique<UnifiedVector>(ndof);

    auto x0 = make_unique<UnifiedVector>(ndof);
    *x0 = sol;

    // Device-resident scalars
    auto s_rho       = r->CreateScalar();
    auto s_rhoLast   = r->CreateScalar();
    auto s_vtrstar   = r->CreateScalar();
    auto s_alpha     = r->CreateScalar();
    auto s_neg_alpha = r->CreateScalar();
    auto s_wnorm_sq  = r->CreateScalar();
    auto s_theta     = r->CreateScalar();
    auto s_c         = r->CreateScalar();
    auto s_tau       = r->CreateScalar();
    auto s_tau_sq    = r->CreateScalar();  // tau^2 for WHILE convergence check
    auto s_eta       = r->CreateScalar();
    auto s_coeff     = r->CreateScalar();
    auto s_beta      = r->CreateScalar();
    auto s_beta_sq   = r->CreateScalar();

    auto& urho       = dynamic_cast<UnifiedScalar&>(*s_rho);
    auto& urhoLast   = dynamic_cast<UnifiedScalar&>(*s_rhoLast);
    auto& uvtrstar   = dynamic_cast<UnifiedScalar&>(*s_vtrstar);
    auto& ualpha     = dynamic_cast<UnifiedScalar&>(*s_alpha);
    auto& uneg_alpha = dynamic_cast<UnifiedScalar&>(*s_neg_alpha);
    auto& uwnorm_sq  = dynamic_cast<UnifiedScalar&>(*s_wnorm_sq);
    auto& utheta     = dynamic_cast<UnifiedScalar&>(*s_theta);
    auto& uc         = dynamic_cast<UnifiedScalar&>(*s_c);
    auto& utau       = dynamic_cast<UnifiedScalar&>(*s_tau);
    auto& utau_sq    = dynamic_cast<UnifiedScalar&>(*s_tau_sq);
    auto& ueta       = dynamic_cast<UnifiedScalar&>(*s_eta);
    auto& ucoeff     = dynamic_cast<UnifiedScalar&>(*s_coeff);
    auto& ubeta      = dynamic_cast<UnifiedScalar&>(*s_beta);
    auto& ubeta_sq   = dynamic_cast<UnifiedScalar&>(*s_beta_sq);

    // Reset vectors and scalars to initial state
    auto InitState = [&]() {
        *x = *x0;
        *r = rhs;
        ga->MultAdd(-1.0, *x0, *r);   // r = rhs - A*x0
        *u = *r;
        *w = *r;
        *d = 0.0;
        ga->Mult(*r, *tmp);
        gc->Mult(*tmp, *v);
        *uhat = *v;
        (*r).InnerProduct(*r, *s_rho);
        SyncNGSStream();
        ualpha     = Div(Const(1.0), Const(1.0));    // 1.0
        uneg_alpha = Neg(Const(1.0));                 // -1.0
        utheta     = Prod(Const(0.0), Const(1.0));   // 0.0
        ueta       = Prod(Const(0.0), Const(1.0));   // 0.0
        ucoeff     = Prod(Const(0.0), Const(1.0));   // 0.0
        utau       = Sqrt(Scal(urho));               // r0norm
        utau_sq    = Scal(urho);                     // r0norm^2 = rho
        urhoLast   = Scal(urho);
        SyncNGSStream();
    };

    InitState();
    double r0norm = utau.GetD();
    if (r0norm == 0.0) { return; }
    double tol = GetPrecision() * r0norm;

    // Warm-up
    ga->Mult(*u, *tmp);
    SyncNGSStream();

    r->UpdateDevice(); u->UpdateDevice(); v->UpdateDevice();
    w->UpdateDevice(); uhat->UpdateDevice(); d->UpdateDevice();
    x->UpdateDevice(); tmp->UpdateDevice();
    SyncNGSStream();

    // Even half-step body. update_tau_sq unused — batch kernels always write tau_sq.
    auto CaptureEven = [&](bool) {
        (*r).InnerProduct(*v, *s_vtrstar);
        ngs_cuda::TFQMREvenBatch1(urho.DevPtr(), uvtrstar.DevPtr(),
                                   utheta.DevPtr(), ueta.DevPtr(),
                                   ualpha.DevPtr(), uneg_alpha.DevPtr(), ucoeff.DevPtr());
        w->Add(*s_neg_alpha, *uhat);
        d->Scale(*s_coeff);
        d->Add(1.0, *u);
        (*w).InnerProduct(*w, *s_wnorm_sq);
        ngs_cuda::TFQMREvenTauBatch(uwnorm_sq.DevPtr(), utau.DevPtr(), ualpha.DevPtr(),
                                     urho.DevPtr(), utheta.DevPtr(), uc.DevPtr(),
                                     utau.DevPtr(), utau_sq.DevPtr(),
                                     ueta.DevPtr(), urhoLast.DevPtr());
        gc->Mult(*d, *tmp);
        x->Add(*s_eta, *tmp);
        u->Add(*s_neg_alpha, *v);  // u = uNext  (d used old u above)
        ga->Mult(*u, *tmp);
        gc->Mult(*tmp, *uhat);
    };

    auto CaptureOdd = [&](bool) {
        ngs_cuda::TFQMROddCoeff(utheta.DevPtr(), ueta.DevPtr(), ualpha.DevPtr(), ucoeff.DevPtr());
        w->Add(*s_neg_alpha, *uhat);
        d->Scale(*s_coeff);
        d->Add(1.0, *u);
        (*w).InnerProduct(*w, *s_wnorm_sq);
        ngs_cuda::TFQMROddTauBatch(uwnorm_sq.DevPtr(), utau.DevPtr(), ualpha.DevPtr(),
                                    utheta.DevPtr(), uc.DevPtr(),
                                    utau.DevPtr(), utau_sq.DevPtr(), ueta.DevPtr());
        gc->Mult(*d, *tmp);
        x->Add(*s_eta, *tmp);
        (*w).InnerProduct(*r, *s_rho);
        ngs_cuda::TFQMROddBeta(urho.DevPtr(), urhoLast.DevPtr(), ubeta.DevPtr(), ubeta_sq.DevPtr());
        u->Scale(*s_beta);
        u->Add(1.0, *w);
        v->Scale(*s_beta_sq);
        v->Add(*s_beta, *uhat);
        ga->Mult(*u, *tmp);
        gc->Mult(*tmp, *uhat);
        v->Add(1.0, *uhat);
    };

    // Graph capture: try WHILE first, fall back to per-iteration 
    ngs_cuda::CudaWhileGraph g_while;
    CudaGraph g_pair;                // even+odd combined — body for WHILE
    CudaGraph g_even, g_odd;         // alternating per-iteration fallback
    bool use_while_graph = false;
    bool use_graph       = false;
    int* iter_count_dev  = nullptr;

    if (!getenv("NO_CUDA_GRAPH")) {
        // Try WHILE graph with single launch for entire solve
        try {
            g_pair.BeginCapture();
            CaptureEven(true);   // update tau_sq in even body
            CaptureOdd(true);    // update tau_sq in odd body
            g_pair.EndCapture();
            SyncNGSStream();

            cudaMalloc(&iter_count_dev, sizeof(int));
            cudaMemset(iter_count_dev, 0, sizeof(int));
            // Each g_pair = 2 TFQMR half-steps; convergence kernel checks sqrt(tau_sq) <= tol
            g_while.Build(g_pair.GetGraph(), utau_sq.DevPtr(), tol,
                          iter_count_dev, (GetMaxSteps() + 1) / 2);
            SyncNGSStream();
            use_while_graph = g_while.IsValid();
        } catch (ngstd::Exception& e) {
            cerr << "[DevTFQMRSolver] WHILE graph capture failed: " << e.What() << endl;
        } catch (...) {
            cerr << "[DevTFQMRSolver] WHILE graph capture failed (unknown exception)" << endl;
        }

        // Fall back to per-iteration even/odd graphs
        if (!use_while_graph) {
            try {
                g_even.BeginCapture();
                CaptureEven(false);
                g_even.EndCapture();
                SyncNGSStream();

                g_odd.BeginCapture();
                CaptureOdd(false);
                g_odd.EndCapture();
                SyncNGSStream();

                use_graph = g_even.IsValid() && g_odd.IsValid();
            } catch (ngstd::Exception& e) {
                cerr << "[DevTFQMRSolver] per-iter graph capture failed: " << e.What() << endl;
            } catch (...) {
                cerr << "[DevTFQMRSolver] per-iter graph capture failed (unknown exception)" << endl;
            }
        }
    }

    // Re-initialise after warm-up / capture
    InitState();

    // Execute 
    if (use_while_graph) {
        cudaMemset(iter_count_dev, 0, sizeof(int));
        g_while.Launch();
        SyncNGSStream();
        cublasSetPointerMode(Get_CuBlas_Handle(), CUBLAS_POINTER_MODE_HOST);
        int pair_count = 0;
        cudaMemcpy(&pair_count, iter_count_dev, sizeof(int), cudaMemcpyDeviceToHost);
        cudaFree(iter_count_dev);
        steps = pair_count * 2;
        sol = *x;
        return;
    }

    // Main loop (per-iteration graph or no-graph)
    for (int iter = 0; iter < GetMaxSteps(); iter++) {
        bool even = (iter % 2 == 0);

        if (use_graph) {
            if (even) g_even.Launch(); else g_odd.Launch();
            SyncNGSStream();
        } else {
            if (even)
                CaptureEven(false);
            else
                CaptureOdd(false);
            SyncNGSStream();
        }

        double tau = utau.GetD();
        if (printrates)
            cout << "TFQMR iter " << iter << "  tau=" << tau << endl;
        if (tau < tol) {
            sol = *x;
            if (printrates)
                cout << "TFQMR converged after " << iter+1 << " iters" << endl;
            return;
        }
    }

    sol = *x;
}

} // namespace ngla
