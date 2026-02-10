/**
 * @file cg_solver.hpp
 * @brief Conjugate Gradient (CG) solver
 */

#ifndef SPARSESOLV_SOLVERS_CG_SOLVER_HPP
#define SPARSESOLV_SOLVERS_CG_SOLVER_HPP

#include "iterative_solver.hpp"

namespace sparsesolv {

/**
 * @brief Preconditioned Conjugate Gradient (PCG) solver
 *
 * Solves symmetric positive definite linear systems Ax = b using the
 * preconditioned conjugate gradient method.
 *
 * The algorithm:
 * 1. r_0 = b - A*x_0
 * 2. z_0 = M^{-1} * r_0
 * 3. p_0 = z_0
 * 4. For k = 0, 1, 2, ...:
 *    a. alpha_k = (r_k, z_k) / (p_k, A*p_k)
 *    b. x_{k+1} = x_k + alpha_k * p_k
 *    c. r_{k+1} = r_k - alpha_k * A*p_k
 *    d. Check convergence
 *    e. z_{k+1} = M^{-1} * r_{k+1}
 *    f. beta_k = (r_{k+1}, z_{k+1}) / (r_k, z_k)
 *    g. p_{k+1} = z_{k+1} + beta_k * p_k
 *
 * Requirements:
 * - Matrix A must be symmetric positive definite
 * - Preconditioner M should also be SPD
 *
 * Example:
 * @code
 * CGSolver<double> solver;
 * solver.config().tolerance = 1e-10;
 * solver.config().max_iterations = 1000;
 *
 * ICPreconditioner<double> precond(1.05);
 * precond.setup(A);
 *
 * auto result = solver.solve(A, b, x, n, &precond);
 * if (result.converged) {
 *     std::cout << "Converged in " << result.iterations << " iterations\n";
 * }
 * @endcode
 *
 * @tparam Scalar The scalar type (double or complex<double>)
 */
template<typename Scalar = double>
class CGSolver : public IterativeSolver<Scalar> {
public:
    std::string name() const override { return "CG"; }

protected:
    SolverResult do_iterate() override {
        const index_t n = this->size_;
        auto& r = this->r_;
        auto& z = this->z_;
        auto& p = this->p_;
        auto& Ap = this->Ap_;
        Scalar* x = this->x_;
        const auto& config = this->config_;

        // Apply preconditioner: z = M^{-1} * r
        this->apply_preconditioner();

        // p = z
        std::copy(z.begin(), z.end(), p.begin());

        // rz_old = (r, z)
        Scalar rz_old = this->dot_product(r.data(), z.data(), n);

        // Main iteration loop
        for (int iter = 0; iter < config.max_iterations; ++iter) {
            // Ap = A * p
            this->A_->multiply(p.data(), Ap.data());

            // alpha = (r, z) / (p, Ap)
            Scalar pAp = this->dot_product(p.data(), Ap.data(), n);

            // Avoid division by zero
            if (std::abs(pAp) < 1e-30) {
                // Numerical breakdown - check if already converged
                double norm_r = this->compute_norm(r.data(), n);
                if (this->check_convergence(norm_r, iter)) {
                    return this->build_result(true, iter + 1, norm_r);
                }
                return this->build_result(false, iter, norm_r);
            }

            Scalar alpha = rz_old / pAp;

            // x = x + alpha * p
            // r = r - alpha * Ap
            #pragma omp parallel for
            for (index_t i = 0; i < n; ++i) {
                x[i] += alpha * p[i];
                r[i] -= alpha * Ap[i];
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

            // z = M^{-1} * r
            this->apply_preconditioner();

            // beta = (r_new, z_new) / (r_old, z_old)
            Scalar rz_new = this->dot_product(r.data(), z.data(), n);
            Scalar beta = rz_new / rz_old;
            rz_old = rz_new;

            // p = z + beta * p
            #pragma omp parallel for
            for (index_t i = 0; i < n; ++i) {
                p[i] = z[i] + beta * p[i];
            }
        }

        // Max iterations reached
        double final_norm = this->compute_norm(r.data(), n);
        return this->build_result(false, config.max_iterations, final_norm);
    }
};

// Type aliases
using CGSolverD = CGSolver<double>;
using CGSolverC = CGSolver<complex_t>;

} // namespace sparsesolv

#endif // SPARSESOLV_SOLVERS_CG_SOLVER_HPP
