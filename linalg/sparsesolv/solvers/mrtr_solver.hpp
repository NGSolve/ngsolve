/**
 * @file mrtr_solver.hpp
 * @brief MRTR (Modified Residual-based Tri-diagonal Reduction) solver
 */

#ifndef SPARSESOLV_SOLVERS_MRTR_SOLVER_HPP
#define SPARSESOLV_SOLVERS_MRTR_SOLVER_HPP

#include "iterative_solver.hpp"

namespace sparsesolv {

/**
 * @brief MRTR (Modified Residual-based Tri-diagonal Reduction) solver
 *
 * MRTR is a Krylov subspace method that is mathematically equivalent to
 * the conjugate residual (CR) method but uses a different recurrence
 * formula. It is suitable for symmetric positive definite systems and
 * can be more stable than CG for some problems.
 *
 * The algorithm (Tsuburaya et al.):
 * 1. r_0 = b - A*x_0
 * 2. u_0 = M^{-1} * r_0
 * 3. y_0 = -r_0, z_0 = M^{-1} * r_0
 * 4. For k = 0, 1, 2, ...:
 *    a. v = A * u
 *    b. w = M^{-1} * v
 *    c. Compute zeta_k, eta_k based on inner products
 *    d. Update p, x, y, r, z, u
 *    e. Check convergence
 *
 * Features:
 * - Minimizes residual norm at each step (unlike CG which minimizes energy)
 * - Works with IC preconditioning (ICMRTR)
 * - Also works with SGS preconditioning (SGS-MRTR)
 *
 * Reference:
 * - Tsuburaya et al., "Convergence acceleration of iterative solvers
 *   for finite element electromagnetic analysis"
 *
 * @tparam Scalar The scalar type (double or complex<double>)
 */
template<typename Scalar = double>
class MRTRSolver : public IterativeSolver<Scalar> {
public:
    std::string name() const override { return "MRTR"; }

protected:
    void allocate_work_vectors() override {
        // Call base class allocation
        IterativeSolver<Scalar>::allocate_work_vectors();

        // Allocate MRTR-specific vectors
        u_.resize(this->size_);
        y_.resize(this->size_);
        v_.resize(this->size_);
        w_.resize(this->size_);
    }

    SolverResult do_iterate() override {
        const index_t n = this->size_;
        auto& r = this->r_;
        auto& z = this->z_;
        auto& p = this->p_;
        Scalar* x = this->x_;
        const auto& config = this->config_;

        // Initialize p to zero
        std::fill(p.begin(), p.end(), Scalar(0));

        // u_0 = M^{-1} * r_0
        this->apply_preconditioner();  // z = M^{-1} * r
        std::copy(z.begin(), z.end(), u_.begin());

        // y_0 = -r_0
        for (index_t i = 0; i < n; ++i) {
            y_[i] = -r[i];
        }

        // z_0 = M^{-1} * r_0 (already computed as z above)

        Scalar zeta = Scalar(1);
        Scalar zeta_old = Scalar(1);
        Scalar eta = Scalar(0);
        Scalar nu = Scalar(1);

        // Main iteration loop
        for (int iter = 0; iter < config.max_iterations; ++iter) {
            // v = A * u
            this->A_->multiply(u_.data(), v_.data());

            // w = M^{-1} * v
            if (this->precond_) {
                this->precond_->apply(v_.data(), w_.data(), n);
            } else {
                std::copy(v_.begin(), v_.end(), w_.begin());
            }

            // Compute inner products
            Scalar w_r = this->dot_product(w_.data(), r.data(), n);
            Scalar v_w = this->dot_product(v_.data(), w_.data(), n);

            if (iter == 0) {
                // First iteration: simplified formulas
                // zeta_0 = (w, r) / (v, w)
                if (std::abs(v_w) < 1e-30) {
                    return this->build_result(false, iter, this->compute_norm(r.data(), n));
                }
                zeta = w_r / v_w;
                zeta_old = zeta;
                eta = Scalar(0);
            } else {
                // General iteration
                Scalar w_y = this->dot_product(w_.data(), y_.data(), n);

                // Denominator for zeta and eta
                Scalar denom = nu * v_w - w_y * w_y;
                if (std::abs(denom) < 1e-30) {
                    return this->build_result(false, iter, this->compute_norm(r.data(), n));
                }

                Scalar inv_denom = Scalar(1) / denom;

                // zeta_k = nu * (w, r) / denom
                zeta = nu * w_r * inv_denom;

                // eta_k = -(w, y) * (w, r) / denom
                eta = -w_y * w_r * inv_denom;
            }

            // nu_{k+1} = zeta * (w, r)
            nu = zeta * w_r;

            // p_k = u + (eta * zeta_old / zeta) * p
            Scalar coeff = eta * zeta_old / zeta;
            #pragma omp parallel for
            for (index_t i = 0; i < n; ++i) {
                p[i] = u_[i] + coeff * p[i];
            }
            zeta_old = zeta;

            // x_{k+1} = x_k + zeta * p
            #pragma omp parallel for
            for (index_t i = 0; i < n; ++i) {
                x[i] += zeta * p[i];
            }

            // y_{k+1} = eta * y + zeta * v
            #pragma omp parallel for
            for (index_t i = 0; i < n; ++i) {
                y_[i] = eta * y_[i] + zeta * v_[i];
            }

            // r_{k+1} = r_k - y_{k+1}
            #pragma omp parallel for
            for (index_t i = 0; i < n; ++i) {
                r[i] -= y_[i];
            }

            // Compute residual norm
            double norm_r = this->compute_norm(r.data(), n);

            // Check convergence
            if (this->check_convergence(norm_r, iter)) {
                return this->build_result(true, iter + 1, norm_r);
            }

            // Check divergence
            if (this->is_diverging()) {
                return this->build_result(false, iter + 1, norm_r);
            }

            // z_{k+1} = eta * z + zeta * w
            #pragma omp parallel for
            for (index_t i = 0; i < n; ++i) {
                z[i] = eta * z[i] + zeta * w_[i];
            }

            // u_{k+1} = u_k - z_{k+1}
            #pragma omp parallel for
            for (index_t i = 0; i < n; ++i) {
                u_[i] -= z[i];
            }
        }

        // Max iterations reached
        double final_norm = this->compute_norm(r.data(), n);
        return this->build_result(false, config.max_iterations, final_norm);
    }

private:
    // MRTR-specific work vectors
    std::vector<Scalar> u_;  // Preconditioned direction
    std::vector<Scalar> y_;  // Residual update direction
    std::vector<Scalar> v_;  // A * u
    std::vector<Scalar> w_;  // M^{-1} * v
};

// Type aliases
using MRTRSolverD = MRTRSolver<double>;
using MRTRSolverC = MRTRSolver<complex_t>;

} // namespace sparsesolv

#endif // SPARSESOLV_SOLVERS_MRTR_SOLVER_HPP
