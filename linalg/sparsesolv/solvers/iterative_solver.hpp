/**
 * @file iterative_solver.hpp
 * @brief Base class for iterative linear solvers
 */

#ifndef SPARSESOLV_SOLVERS_ITERATIVE_SOLVER_HPP
#define SPARSESOLV_SOLVERS_ITERATIVE_SOLVER_HPP

#include "../core/types.hpp"
#include "../core/solver_config.hpp"
#include "../core/sparse_matrix_view.hpp"
#include "../core/preconditioner.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace sparsesolv {

/**
 * @brief Abstract base class for iterative linear solvers
 *
 * This class provides a common interface for iterative methods such as CG,
 * MRTR, BiCGSTAB, etc. It handles:
 * - Configuration management
 * - Convergence checking
 * - Divergence detection
 * - Residual history tracking
 * - Best result saving
 *
 * To implement a new solver:
 * 1. Override do_iterate() to implement the specific iteration scheme
 * 2. Optionally override name() to return a descriptive name
 *
 * @tparam Scalar The scalar type (double or complex<double>)
 */
template<typename Scalar = double>
class IterativeSolver {
public:
    using value_type = Scalar;

    virtual ~IterativeSolver() = default;

    /**
     * @brief Solve the linear system Ax = b
     *
     * @param A The system matrix
     * @param b Right-hand side vector
     * @param x Solution vector (initial guess on input, solution on output)
     * @param size System size
     * @param precond Preconditioner (optional, nullptr for no preconditioning)
     * @return SolverResult containing convergence info
     */
    SolverResult solve(
        const SparseMatrixView<Scalar>& A,
        const Scalar* b,
        Scalar* x,
        index_t size,
        const Preconditioner<Scalar>* precond = nullptr
    ) {
        // Initialize
        size_ = size;
        A_ = &A;
        b_ = b;
        x_ = x;
        precond_ = precond;

        // Allocate work vectors
        allocate_work_vectors();

        // Compute initial state
        prepare_iteration();

        // Run the iteration
        return do_iterate();
    }

    /**
     * @brief Solve with std::vector interface
     */
    SolverResult solve(
        const SparseMatrixView<Scalar>& A,
        const std::vector<Scalar>& b,
        std::vector<Scalar>& x,
        const Preconditioner<Scalar>* precond = nullptr
    ) {
        if (x.size() != b.size()) {
            x.resize(b.size());
        }
        return solve(A, b.data(), x.data(), static_cast<index_t>(b.size()), precond);
    }

    /// Set solver configuration
    void set_config(const SolverConfig& config) { config_ = config; }

    /// Get solver configuration
    const SolverConfig& config() const { return config_; }

    /// Get mutable configuration reference
    SolverConfig& config() { return config_; }

    /// Get the name of the solver
    virtual std::string name() const = 0;

protected:
    /**
     * @brief Implement the specific iteration scheme
     *
     * This method should perform the iteration until convergence or
     * max iterations. Use the helper methods for convergence checking.
     *
     * @return SolverResult with final state
     */
    virtual SolverResult do_iterate() = 0;

    /**
     * @brief Allocate work vectors needed by the solver
     *
     * Override this to allocate additional work vectors beyond r_, z_, p_, Ap_
     */
    virtual void allocate_work_vectors() {
        r_.resize(size_);
        z_.resize(size_);
        p_.resize(size_);
        Ap_.resize(size_);
    }

    /**
     * @brief Prepare for iteration (compute initial residual, etc.)
     */
    virtual void prepare_iteration() {
        // Initialize x to zero if requested
        if (config_.save_best_result) {
            best_x_.resize(size_);
            best_residual_ = std::numeric_limits<double>::max();
        }

        // Compute initial residual: r = b - A*x
        A_->multiply(x_, r_.data());
        #pragma omp parallel for
        for (index_t i = 0; i < size_; ++i) {
            r_[i] = b_[i] - r_[i];
        }

        // Compute norm of b for normalization
        norm_b_ = compute_norm(b_, size_);

        // Compute initial residual norm
        double init_norm_r = compute_norm(r_.data(), size_);

        // Set normalizer based on config
        switch (config_.norm_type) {
            case NormType::RHS:
                normalizer_ = norm_b_;
                break;
            case NormType::InitialResidual:
                normalizer_ = init_norm_r;
                break;
            case NormType::Custom:
                normalizer_ = config_.custom_norm;
                break;
        }

        // Avoid division by zero
        if (normalizer_ < 1e-15) {
            normalizer_ = 1.0;
        }

        // Clear residual history
        if (config_.save_residual_history) {
            residual_history_.clear();
            residual_history_.push_back(init_norm_r / normalizer_);
        }

        // Reset divergence counter
        bad_count_ = 0;
    }

    /**
     * @brief Apply preconditioner: z = M^{-1} * r
     */
    void apply_preconditioner() {
        if (precond_) {
            precond_->apply(r_.data(), z_.data(), size_);
        } else {
            // No preconditioning: z = r
            std::copy(r_.begin(), r_.end(), z_.begin());
        }
    }

    /**
     * @brief Check convergence and update tracking
     * @param norm_r Current residual norm (unnormalized)
     * @param iter Current iteration number
     * @return true if converged
     */
    bool check_convergence(double norm_r, int iter) {
        double rel_residual = norm_r / normalizer_;

        // Save to history
        if (config_.save_residual_history) {
            residual_history_.push_back(rel_residual);
        }

        // Update best result
        if (config_.save_best_result && rel_residual < best_residual_) {
            best_residual_ = rel_residual;
            std::copy(x_, x_ + size_, best_x_.begin());
        }

        // Check convergence
        if (rel_residual < config_.tolerance || norm_r < config_.abs_tolerance) {
            return true;
        }

        // Check divergence
        if (config_.divergence_check == DivergenceCheck::StagnationCount) {
            if (rel_residual >= best_residual_ * config_.divergence_threshold) {
                bad_count_++;
                if (bad_count_ >= config_.divergence_count) {
                    return false; // Will be marked as not converged
                }
            } else {
                bad_count_ = 0;
            }
        }

        return false;
    }

    /**
     * @brief Check if we should stop due to divergence
     */
    bool is_diverging() const {
        return config_.divergence_check == DivergenceCheck::StagnationCount &&
               bad_count_ >= config_.divergence_count;
    }

    /**
     * @brief Build the result structure
     */
    SolverResult build_result(bool converged, int iterations, double final_residual) {
        SolverResult result;
        result.converged = converged;
        result.iterations = iterations;
        result.final_residual = final_residual / normalizer_;

        if (config_.save_residual_history) {
            result.residual_history = std::move(residual_history_);
        }

        // If not converged but we have best result, use it
        if (!converged && config_.save_best_result) {
            std::copy(best_x_.begin(), best_x_.end(), x_);
            result.final_residual = best_residual_;
        }

        return result;
    }

    /**
     * @brief Compute Euclidean norm of a vector
     */
    static double compute_norm(const Scalar* v, index_t size) {
        double sum = 0.0;
        #pragma omp parallel for reduction(+:sum)
        for (index_t i = 0; i < size; ++i) {
            sum += std::norm(v[i]); // std::norm for complex gives |z|^2
        }
        return std::sqrt(sum);
    }

    /**
     * @brief Compute dot product of two vectors
     *
     * For real types: standard dot product
     * For complex types: conjugate of first argument (Hermitian inner product)
     */
    static Scalar dot_product(const Scalar* a, const Scalar* b, index_t size) {
        Scalar sum = Scalar(0);
        if constexpr (std::is_same_v<Scalar, complex_t>) {
            // Complex: MSVC OpenMP doesn't support reduction on complex types
            for (index_t i = 0; i < size; ++i) {
                sum += std::conj(a[i]) * b[i];
            }
        } else {
            // Real: standard dot product with OpenMP reduction
            #pragma omp parallel for reduction(+:sum)
            for (index_t i = 0; i < size; ++i) {
                sum += a[i] * b[i];
            }
        }
        return sum;
    }

    // Configuration
    SolverConfig config_;

    // Problem data (set during solve)
    index_t size_ = 0;
    const SparseMatrixView<Scalar>* A_ = nullptr;
    const Scalar* b_ = nullptr;
    Scalar* x_ = nullptr;
    const Preconditioner<Scalar>* precond_ = nullptr;

    // Work vectors
    std::vector<Scalar> r_;   // Residual
    std::vector<Scalar> z_;   // Preconditioned residual
    std::vector<Scalar> p_;   // Search direction
    std::vector<Scalar> Ap_;  // Matrix-vector product A*p

    // Convergence tracking
    double norm_b_ = 0.0;
    double normalizer_ = 1.0;
    std::vector<double> residual_history_;

    // Best result tracking
    std::vector<Scalar> best_x_;
    double best_residual_ = std::numeric_limits<double>::max();

    // Divergence tracking
    int bad_count_ = 0;
};

// Type aliases
using IterativeSolverD = IterativeSolver<double>;
using IterativeSolverC = IterativeSolver<complex_t>;

} // namespace sparsesolv

#endif // SPARSESOLV_SOLVERS_ITERATIVE_SOLVER_HPP
