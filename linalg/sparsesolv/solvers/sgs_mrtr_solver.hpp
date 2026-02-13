/**
 * @file sgs_mrtr_solver.hpp
 * @brief SGS-MRTR solver (MRTR with Symmetric Gauss-Seidel preconditioning)
 *
 * This is a specialized solver that applies the SGS preconditioner using
 * a split formula (L and L^T separately) within the MRTR iteration.
 * This is different from the generic MRTR solver which applies M^{-1} as a whole.
 *
 * Now inherits from IterativeSolver to gain:
 * - save_best_result: tracks and restores best solution during iteration
 * - divergence_check: early termination on stagnation
 * - residual_history: residual tracking per iteration
 */

#ifndef SPARSESOLV_SOLVERS_SGS_MRTR_SOLVER_HPP
#define SPARSESOLV_SOLVERS_SGS_MRTR_SOLVER_HPP

#include "iterative_solver.hpp"
#include "../core/constants.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>

namespace sparsesolv {

/**
 * @brief SGS-MRTR solver
 *
 * Implements the MRTR algorithm with specialized SGS preconditioning using
 * the split formula where L (forward) and L^T (backward) are applied separately.
 *
 * The solver works in a diagonally-scaled space (DAD) where D = diag(A)^{-1/2}.
 * This improves conditioning and is built into the iteration.
 *
 * The algorithm works in the preconditioned residual space:
 * - rd = L^{-1} * r (preconditioned residual)
 * - ARd = L^T * rd + L^{-1} * (rd - L^T * rd) (approximate M^{-1} * A * rd)
 *
 * Inherits from IterativeSolver for save_best_result, divergence_check,
 * and residual_history tracking. The precond parameter in solve() is ignored
 * (SGS preconditioning is built into the iteration).
 *
 * @tparam Scalar The scalar type (double or complex<double>)
 */
template<typename Scalar = double>
class SGSMRTRSolver : public IterativeSolver<Scalar> {
public:
    SGSMRTRSolver() = default;

    std::string name() const override { return "SGS-MRTR"; }

protected:
    void allocate_work_vectors() override {
        const index_t n = this->size_;

        // Base class vectors (r_ used for scaled residual)
        this->r_.resize(n);
        this->z_.resize(n);   // not used but allocated by base
        this->p_.resize(n);
        this->Ap_.resize(n);  // not used but allocated by base

        // SGS-specific vectors
        x2_.resize(n);
        D_.resize(n);
        inv_D_.resize(n);
        b2_.resize(n);
        rd_.resize(n);
        u_.resize(n);
        y_.resize(n);
        Ard_.resize(n);
        temp_.resize(n);
    }

    void prepare_iteration() override {
        const index_t n = this->size_;

        // Save original x pointer (x_ will be redirected to x2_)
        x_orig_ = this->x_;

        // Extract diagonal scaling: D[i] = 1/sqrt(|A[i,i]|)
        #pragma omp parallel for
        for (index_t i = 0; i < n; ++i) {
            Scalar aii = this->A_->diagonal(i);
            if (std::abs(aii) > constants::MIN_DIAGONAL_TOLERANCE) {
                D_[i] = Scalar(1) / std::sqrt(std::abs(aii));
                inv_D_[i] = std::sqrt(std::abs(aii));
            } else {
                D_[i] = Scalar(1);
                inv_D_[i] = Scalar(1);
            }
        }

        // Scaled RHS: b2 = D * b
        #pragma omp parallel for
        for (index_t i = 0; i < n; ++i)
            b2_[i] = D_[i] * this->b_[i];

        // Scaled initial guess: x2 = inv_D * x
        #pragma omp parallel for
        for (index_t i = 0; i < n; ++i)
            x2_[i] = inv_D_[i] * x_orig_[i];

        // Build L (lower tri of DAD) and L^T
        build_L_and_Lt();

        // Compute initial residual in scaled space: r = D*(b - A*x)
        this->A_->multiply(x_orig_, temp_.data());
        #pragma omp parallel for
        for (index_t i = 0; i < n; ++i)
            this->r_[i] = D_[i] * (this->b_[i] - temp_[i]);

        // Normalizer: ||L * b2||
        multiply_L(b2_.data(), temp_.data());
        double norm_Lb = this->compute_norm(temp_.data(), n);
        this->normalizer_ = (norm_Lb > 1e-15) ? norm_Lb : 1.0;

        // Redirect x_ to x2_ so base class best-result tracking works on x2
        this->x_ = x2_.data();

        // Initialize best result tracking
        if (this->config_.save_best_result) {
            this->best_x_.resize(n);
            this->best_residual_ = std::numeric_limits<double>::max();
        }
        this->bad_count_ = 0;

        // Initial residual history
        double norm_r = this->compute_norm(this->r_.data(), n);
        if (this->config_.save_residual_history) {
            this->residual_history_.clear();
            this->residual_history_.push_back(norm_r / this->normalizer_);
        }

        // rd = L^{-1} * r (forward solve)
        forward_solve(this->r_.data(), rd_.data());

        // y_0 = -rd
        #pragma omp parallel for
        for (index_t i = 0; i < n; ++i)
            y_[i] = -rd_[i];

        // Initialize p to zero
        std::fill(this->p_.begin(), this->p_.end(), Scalar(0));
    }

    SolverResult do_iterate() override {
        const index_t n = this->size_;

        // Check if already converged
        double norm_b2 = this->compute_norm(b2_.data(), n);
        if (norm_b2 < constants::BREAKDOWN_THRESHOLD) norm_b2 = 1.0;
        double norm_r = this->compute_norm(this->r_.data(), n);
        if (norm_r / norm_b2 < this->config_.tolerance * 0.1) {
            auto result = this->build_result(true, 0, norm_r);
            convert_x2_to_x();
            return result;
        }

        Scalar zeta(1), zeta_old(1), eta(0), nu(1);

        int iter = 0;
        for (iter = 0; iter < this->config_.max_iterations; ++iter) {
            // u = L^{-T} * rd (backward solve with L^T)
            backward_solve(rd_.data(), u_.data());

            // ARd = u + L^{-1} * (rd - u)
            #pragma omp parallel for
            for (index_t i = 0; i < n; ++i)
                temp_[i] = rd_[i] - u_[i];
            forward_solve(temp_.data(), Ard_.data());
            #pragma omp parallel for
            for (index_t i = 0; i < n; ++i)
                Ard_[i] += u_[i];

            // Compute inner products
            Scalar Ar_r = this->dot_product(Ard_.data(), rd_.data(), n);
            Scalar Ar_Ar = this->dot_product(Ard_.data(), Ard_.data(), n);

            if (iter == 0) {
                zeta = Ar_r / Ar_Ar;
                zeta_old = zeta;
                eta = Scalar(0);
            } else {
                Scalar Ar_y = this->dot_product(Ard_.data(), y_.data(), n);
                Scalar denom = nu * Ar_Ar - Ar_y * Ar_y;

                if (std::abs(denom) < constants::DENOMINATOR_BREAKDOWN) {
                    denom = (std::real(denom) >= 0)
                        ? Scalar(constants::DENOMINATOR_BREAKDOWN)
                        : Scalar(-constants::DENOMINATOR_BREAKDOWN);
                }

                Scalar inv_denom = Scalar(1) / denom;
                zeta = nu * Ar_r * inv_denom;
                eta = -Ar_y * Ar_r * inv_denom;
            }

            nu = zeta * Ar_r;

            // p = u + (eta * zeta_old / zeta) * p
            Scalar coeff = eta * zeta_old / zeta;
            #pragma omp parallel for
            for (index_t i = 0; i < n; ++i)
                this->p_[i] = u_[i] + coeff * this->p_[i];
            zeta_old = zeta;

            // x2 = x2 + zeta * p (x_ points to x2_)
            #pragma omp parallel for
            for (index_t i = 0; i < n; ++i)
                this->x_[i] += zeta * this->p_[i];

            // y = eta * y + zeta * ARd
            #pragma omp parallel for
            for (index_t i = 0; i < n; ++i)
                y_[i] = eta * y_[i] + zeta * Ard_[i];

            // rd = rd - y
            #pragma omp parallel for
            for (index_t i = 0; i < n; ++i)
                rd_[i] -= y_[i];

            // Convergence check (base class handles history, best-result, divergence)
            norm_r = this->compute_norm(rd_.data(), n);
            if (this->check_convergence(norm_r, iter)) {
                auto result = this->build_result(true, iter + 1, norm_r);
                convert_x2_to_x();
                return result;
            }

            if (this->is_diverging()) {
                auto result = this->build_result(false, iter + 1, norm_r);
                convert_x2_to_x();
                return result;
            }
        }

        norm_r = this->compute_norm(rd_.data(), n);
        auto result = this->build_result(false, iter, norm_r);
        convert_x2_to_x();
        return result;
    }

private:
    // Original x pointer (before redirection to x2_)
    Scalar* x_orig_ = nullptr;

    // SGS-specific work vectors
    std::vector<Scalar> x2_;      // Scaled solution: x2 = inv_D * x
    std::vector<Scalar> D_;       // D[i] = 1/sqrt(|A[i,i]|)
    std::vector<Scalar> inv_D_;   // inv_D[i] = sqrt(|A[i,i]|)
    std::vector<Scalar> b2_;      // Scaled RHS: b2 = D * b
    std::vector<Scalar> rd_;      // Preconditioned residual: rd = L^{-1} * r
    std::vector<Scalar> u_;       // Work vector for L^{-T} * rd
    std::vector<Scalar> y_;       // Accumulated correction
    std::vector<Scalar> Ard_;     // Approximate M^{-1} * A * rd
    std::vector<Scalar> temp_;    // Temporary buffer

    // L (lower triangular of DAD) stored as CSR
    std::vector<index_t> L_row_ptr_;
    std::vector<index_t> L_col_idx_;
    std::vector<Scalar> L_values_;

    // L^T stored as CSR
    std::vector<index_t> Lt_row_ptr_;
    std::vector<index_t> Lt_col_idx_;
    std::vector<Scalar> Lt_values_;

    /**
     * @brief Convert scaled solution x2_ back to original space x
     *
     * Must be called AFTER build_result(), because build_result() may
     * restore best_x_ to x_ (= x2_) when not converged.
     * After that, we convert x2_ → x_orig_ via x = D * x2.
     */
    void convert_x2_to_x() {
        #pragma omp parallel for
        for (index_t i = 0; i < this->size_; ++i)
            x_orig_[i] = D_[i] * x2_[i];
    }

    /**
     * @brief Extract lower triangular part of DAD and compute L^T
     */
    void build_L_and_Lt() {
        const index_t n = this->size_;
        const auto& A = *this->A_;

        L_row_ptr_.resize(n + 1);
        L_col_idx_.clear();
        L_values_.clear();

        L_row_ptr_[0] = 0;
        for (index_t i = 0; i < n; ++i) {
            auto [start, end] = A.row_range(i);
            for (index_t k = start; k < end; ++k) {
                index_t j = A.col_idx()[k];
                if (j <= i) {
                    L_col_idx_.push_back(j);
                    L_values_.push_back(D_[i] * A.values()[k] * D_[j]);
                }
            }
            L_row_ptr_[i + 1] = static_cast<index_t>(L_col_idx_.size());
        }

        // Compute L^T (transpose of L)
        const size_t nnz = L_col_idx_.size();
        Lt_row_ptr_.assign(n + 1, 0);
        Lt_col_idx_.resize(nnz);
        Lt_values_.resize(nnz);

        for (size_t k = 0; k < nnz; ++k)
            Lt_row_ptr_[L_col_idx_[k] + 1]++;
        for (index_t i = 0; i < n; ++i)
            Lt_row_ptr_[i + 1] += Lt_row_ptr_[i];

        std::vector<index_t> counter(n, 0);
        for (index_t i = 0; i < n; ++i) {
            for (index_t k = L_row_ptr_[i]; k < L_row_ptr_[i + 1]; ++k) {
                index_t j = L_col_idx_[k];
                index_t pos = Lt_row_ptr_[j] + counter[j];
                Lt_col_idx_[pos] = i;
                Lt_values_[pos] = L_values_[k];
                counter[j]++;
            }
        }
    }

    /**
     * @brief Forward solve: L * y = rhs  →  y = L^{-1} * rhs
     */
    void forward_solve(const Scalar* rhs, Scalar* y) const {
        for (index_t i = 0; i < this->size_; ++i) {
            Scalar s = rhs[i];
            const index_t row_end = L_row_ptr_[i + 1] - 1; // Exclude diagonal
            for (index_t k = L_row_ptr_[i]; k < row_end; ++k)
                s -= L_values_[k] * y[L_col_idx_[k]];
            y[i] = s / L_values_[row_end]; // Divide by diagonal
        }
    }

    /**
     * @brief Backward solve: L^T * y = rhs  →  y = L^{-T} * rhs
     */
    void backward_solve(const Scalar* rhs, Scalar* y) const {
        for (index_t i = this->size_; i-- > 0;) {
            Scalar s = rhs[i];
            const index_t row_end = Lt_row_ptr_[i + 1];
            for (index_t k = Lt_row_ptr_[i] + 1; k < row_end; ++k)
                s -= Lt_values_[k] * y[Lt_col_idx_[k]];
            y[i] = s / Lt_values_[Lt_row_ptr_[i]]; // Divide by diagonal
        }
    }

    /**
     * @brief Matrix-vector multiplication: y = L * x
     */
    void multiply_L(const Scalar* x, Scalar* y) const {
        #pragma omp parallel for
        for (index_t i = 0; i < this->size_; ++i) {
            Scalar s = Scalar(0);
            for (index_t k = L_row_ptr_[i]; k < L_row_ptr_[i + 1]; ++k)
                s += L_values_[k] * x[L_col_idx_[k]];
            y[i] = s;
        }
    }
};

// Type aliases
using SGSMRTRSolverD = SGSMRTRSolver<double>;
using SGSMRTRSolverC = SGSMRTRSolver<complex_t>;

} // namespace sparsesolv

#endif // SPARSESOLV_SOLVERS_SGS_MRTR_SOLVER_HPP
