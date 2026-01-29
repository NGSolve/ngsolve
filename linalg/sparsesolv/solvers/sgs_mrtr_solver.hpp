/**
 * @file sgs_mrtr_solver.hpp
 * @brief SGS-MRTR solver (MRTR with Symmetric Gauss-Seidel preconditioning)
 *
 * This is a specialized solver that applies the SGS preconditioner using
 * a split formula (L and L^T separately) within the MRTR iteration.
 * This is different from the generic MRTR solver which applies M^{-1} as a whole.
 */

#ifndef SPARSESOLV_SOLVERS_SGS_MRTR_SOLVER_HPP
#define SPARSESOLV_SOLVERS_SGS_MRTR_SOLVER_HPP

#include "../core/sparse_matrix_view.hpp"
#include "../core/solver_config.hpp"
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
 * The algorithm works in the preconditioned residual space:
 * - rd = L^{-1} * r (preconditioned residual)
 * - ARd = L^T * rd + L^{-1} * (rd - L^T * rd) (approximate M^{-1} * A * rd)
 *
 * @tparam Scalar The scalar type (double or complex<double>)
 */
template<typename Scalar = double>
class SGSMRTRSolver {
public:
    SGSMRTRSolver() = default;

    /// Set solver configuration
    void set_config(const SolverConfig& config) { config_ = config; }

    /// Get solver configuration (mutable)
    SolverConfig& config() { return config_; }

    /// Get solver configuration (const)
    const SolverConfig& config() const { return config_; }

    /**
     * @brief Solve Ax = b using SGS-MRTR
     *
     * @param A The system matrix (must be symmetric positive definite)
     * @param b Right-hand side vector
     * @param x Solution vector (initial guess, modified in place)
     * @param size System size
     * @return SolverResult with convergence information
     */
    SolverResult solve(
        const SparseMatrixView<Scalar>& A,
        const Scalar* b,
        Scalar* x,
        index_t size
    ) {
        size_ = size;

        // Extract diagonal and compute D = diag(A)^{-1/2} for scaling
        std::vector<Scalar> D(size);       // D[i] = 1/sqrt(A[i,i])
        std::vector<Scalar> inv_D(size);   // inv_D[i] = sqrt(A[i,i])
        for (index_t i = 0; i < size; ++i) {
            Scalar aii = A.diagonal(i);
            if (std::abs(aii) > 1e-15) {
                D[i] = Scalar(1) / std::sqrt(std::abs(aii));
                inv_D[i] = std::sqrt(std::abs(aii));
            } else {
                D[i] = Scalar(1);
                inv_D[i] = Scalar(1);
            }
        }

        // Compute scaled RHS: b2 = D * b
        std::vector<Scalar> b2(size);
        for (index_t i = 0; i < size; ++i) {
            b2[i] = D[i] * b[i];
        }

        // Compute scaled initial guess: x2 = D^{-1} * x
        std::vector<Scalar> x2(size);
        for (index_t i = 0; i < size; ++i) {
            x2[i] = inv_D[i] * x[i];
        }

        // Extract lower triangular part of DAD (stored as CSR)
        // DAD[i,j] = D[i] * A[i,j] * D[j]
        std::vector<index_t> L_row_ptr(size + 1);
        std::vector<index_t> L_col_idx;
        std::vector<Scalar> L_values;

        L_row_ptr[0] = 0;
        for (index_t i = 0; i < size; ++i) {
            auto [start, end] = A.row_range(i);
            for (index_t k = start; k < end; ++k) {
                index_t j = A.col_idx()[k];
                if (j <= i) {
                    L_col_idx.push_back(j);
                    // DAD[i,j] = D[i] * A[i,j] * D[j]
                    L_values.push_back(D[i] * A.values()[k] * D[j]);
                }
            }
            L_row_ptr[i + 1] = static_cast<index_t>(L_col_idx.size());
        }

        // Compute L^T (transpose of L)
        std::vector<index_t> Lt_row_ptr(size + 1, 0);
        std::vector<index_t> Lt_col_idx(L_col_idx.size());
        std::vector<Scalar> Lt_values(L_values.size());

        // Count entries per row of L^T
        for (size_t k = 0; k < L_col_idx.size(); ++k) {
            Lt_row_ptr[L_col_idx[k] + 1]++;
        }
        for (index_t i = 0; i < size; ++i) {
            Lt_row_ptr[i + 1] += Lt_row_ptr[i];
        }

        // Fill L^T
        std::vector<index_t> counter(size, 0);
        for (index_t i = 0; i < size; ++i) {
            for (index_t k = L_row_ptr[i]; k < L_row_ptr[i + 1]; ++k) {
                index_t j = L_col_idx[k];
                index_t pos = Lt_row_ptr[j] + counter[j];
                Lt_col_idx[pos] = i;
                Lt_values[pos] = L_values[k];
                counter[j]++;
            }
        }

        // Allocate work vectors
        std::vector<Scalar> r(size);    // residual
        std::vector<Scalar> rd(size);   // preconditioned residual
        std::vector<Scalar> p(size, Scalar(0));
        std::vector<Scalar> u(size);
        std::vector<Scalar> y(size);
        std::vector<Scalar> Ard(size);
        std::vector<Scalar> temp(size);

        // Compute initial residual in scaled space: r = D * (b - A*x)
        // Since DAD * x2 = D * A * D * (D^{-1} * x) = D * A * x
        // We have: r = b2 - DAD * x2 = D * b - D * A * x = D * (b - A*x)
        A.multiply(x, temp.data());  // temp = A*x
        for (index_t i = 0; i < size; ++i) {
            r[i] = D[i] * (b[i] - temp[i]);  // r = D * (b - A*x)
            u[i] = b2[i];
        }

        // Compute norm of b2
        double norm_b = compute_norm(b2.data(), size);
        if (norm_b < 1e-30) {
            norm_b = 1.0;
        }

        // Check if already converged
        double norm_r = compute_norm(r.data(), size);
        if (norm_r / norm_b < config_.tolerance * 0.1) {
            // Convert x2 back to x: x = D * x2
            for (index_t i = 0; i < size; ++i) {
                x[i] = D[i] * x2[i];
            }
            return build_result(true, 0, norm_r);
        }

        // Compute normalized ||L * b2|| for convergence check
        // temp = L * b2 (matrix-vector multiplication)
        multiply_L(L_row_ptr, L_col_idx, L_values, u.data(), temp.data());
        double norm_Lb = compute_norm(temp.data(), size);

        // rd = L^{-1} * r (forward solve)
        forward_solve(L_row_ptr, L_col_idx, L_values, r.data(), rd.data());

        // y_0 = -rd
        for (index_t i = 0; i < size; ++i) {
            y[i] = -rd[i];
        }

        // Store residual history if requested
        if (config_.save_residual_history) {
            residual_history_.clear();
            residual_history_.push_back(norm_r / norm_b);
        }

        Scalar zeta = Scalar(1);
        Scalar zeta_old = Scalar(1);
        Scalar eta = Scalar(0);
        Scalar nu = Scalar(1);

        // Main iteration loop (works in scaled space)
        int iter = 0;
        for (iter = 0; iter < config_.max_iterations; ++iter) {
            // u = L^{-T} * rd (backward solve with L^T)
            backward_solve(Lt_row_ptr, Lt_col_idx, Lt_values, rd.data(), u.data());

            // ARd = u + L^{-1} * (rd - u)
            for (index_t i = 0; i < size; ++i) {
                temp[i] = rd[i] - u[i];
            }
            forward_solve(L_row_ptr, L_col_idx, L_values, temp.data(), Ard.data());
            for (index_t i = 0; i < size; ++i) {
                Ard[i] += u[i];
            }

            // Compute inner products
            Scalar Ar_r = dot_product(Ard.data(), rd.data(), size);
            Scalar Ar_Ar = dot_product(Ard.data(), Ard.data(), size);

            if (iter == 0) {
                // First iteration: simplified formulas
                zeta = Ar_r / Ar_Ar;
                zeta_old = zeta;
                eta = Scalar(0);
            } else {
                // General iteration
                Scalar Ar_y = dot_product(Ard.data(), y.data(), size);

                // Denominator for zeta and eta
                Scalar denom = nu * Ar_Ar - Ar_y * Ar_y;

                // Regularize if denominator is too small (numerical stability)
                if (std::abs(denom) < 1e-60) {
                    denom = (denom >= 0 ? 1e-60 : -1e-60);
                }

                Scalar inv_denom = Scalar(1) / denom;
                zeta = nu * Ar_r * inv_denom;
                eta = -Ar_y * Ar_r * inv_denom;
            }

            // nu_{k+1} = zeta * (ARd, rd)
            nu = zeta * Ar_r;

            // p = u + (eta * zeta_old / zeta) * p
            Scalar coeff = eta * zeta_old / zeta;
            for (index_t i = 0; i < size; ++i) {
                p[i] = u[i] + coeff * p[i];
            }
            zeta_old = zeta;

            // x2 = x2 + zeta * p (update in scaled space)
            for (index_t i = 0; i < size; ++i) {
                x2[i] += zeta * p[i];
            }

            // y = eta * y + zeta * ARd
            for (index_t i = 0; i < size; ++i) {
                y[i] = eta * y[i] + zeta * Ard[i];
            }

            // rd = rd - y
            for (index_t i = 0; i < size; ++i) {
                rd[i] -= y[i];
            }

            // Convergence check
            norm_r = compute_norm(rd.data(), size);
            double relative_residual = norm_r / norm_Lb;

            if (config_.save_residual_history) {
                residual_history_.push_back(relative_residual);
            }

            if (relative_residual < config_.tolerance) {
                // Convert x2 back to x: x = D * x2
                for (index_t i = 0; i < size; ++i) {
                    x[i] = D[i] * x2[i];
                }
                return build_result(true, iter + 1, relative_residual);
            }
        }

        // Convert x2 back to x: x = D * x2
        for (index_t i = 0; i < size; ++i) {
            x[i] = D[i] * x2[i];
        }
        return build_result(false, iter, norm_r / norm_Lb);
    }

    /**
     * @brief Solve with std::vector interface
     */
    SolverResult solve(
        const SparseMatrixView<Scalar>& A,
        const std::vector<Scalar>& b,
        std::vector<Scalar>& x
    ) {
        if (x.size() != b.size()) {
            x.resize(b.size());
        }
        return solve(A, b.data(), x.data(), static_cast<index_t>(b.size()));
    }

    std::string name() const { return "SGS-MRTR"; }

private:
    SolverConfig config_;
    index_t size_ = 0;
    std::vector<double> residual_history_;

    /**
     * @brief Forward solve: L * y = x
     */
    void forward_solve(
        const std::vector<index_t>& row_ptr,
        const std::vector<index_t>& col_idx,
        const std::vector<Scalar>& values,
        const Scalar* x,
        Scalar* y
    ) const {
        for (index_t i = 0; i < size_; ++i) {
            Scalar s = x[i];
            const index_t row_end = row_ptr[i + 1] - 1;  // Exclude diagonal
            for (index_t k = row_ptr[i]; k < row_end; ++k) {
                s -= values[k] * y[col_idx[k]];
            }
            // Divide by diagonal (last element in row)
            y[i] = s / values[row_end];
        }
    }

    /**
     * @brief Backward solve: L^T * y = x
     */
    void backward_solve(
        const std::vector<index_t>& row_ptr,
        const std::vector<index_t>& col_idx,
        const std::vector<Scalar>& values,
        const Scalar* x,
        Scalar* y
    ) const {
        for (index_t i = size_; i-- > 0;) {
            Scalar s = x[i];
            const index_t row_end = row_ptr[i + 1];
            // Skip diagonal (first element)
            for (index_t k = row_ptr[i] + 1; k < row_end; ++k) {
                s -= values[k] * y[col_idx[k]];
            }
            // Divide by diagonal (first element in row)
            y[i] = s / values[row_ptr[i]];
        }
    }

    /**
     * @brief Matrix-vector multiplication: y = L * x
     */
    void multiply_L(
        const std::vector<index_t>& row_ptr,
        const std::vector<index_t>& col_idx,
        const std::vector<Scalar>& values,
        const Scalar* x,
        Scalar* y
    ) const {
        for (index_t i = 0; i < size_; ++i) {
            Scalar s = Scalar(0);
            for (index_t k = row_ptr[i]; k < row_ptr[i + 1]; ++k) {
                s += values[k] * x[col_idx[k]];
            }
            y[i] = s;
        }
    }

    static double compute_norm(const Scalar* v, index_t n) {
        double sum = 0.0;
        for (index_t i = 0; i < n; ++i) {
            sum += std::abs(v[i]) * std::abs(v[i]);
        }
        return std::sqrt(sum);
    }

    static Scalar dot_product(const Scalar* a, const Scalar* b, index_t n) {
        Scalar sum = Scalar(0);
        if constexpr (std::is_same_v<Scalar, complex_t>) {
            for (index_t i = 0; i < n; ++i) {
                sum += std::conj(a[i]) * b[i];
            }
        } else {
            for (index_t i = 0; i < n; ++i) {
                sum += a[i] * b[i];
            }
        }
        return sum;
    }

    SolverResult build_result(bool converged, int iterations, double final_residual) const {
        SolverResult result;
        result.converged = converged;
        result.iterations = iterations;
        result.final_residual = final_residual;
        result.residual_history = residual_history_;
        return result;
    }
};

// Type aliases
using SGSMRTRSolverD = SGSMRTRSolver<double>;
using SGSMRTRSolverC = SGSMRTRSolver<complex_t>;

} // namespace sparsesolv

#endif // SPARSESOLV_SOLVERS_SGS_MRTR_SOLVER_HPP
