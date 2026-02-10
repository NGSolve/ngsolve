/**
 * @file solver_config.hpp
 * @brief Unified configuration for all SparseSolv solvers
 */

#ifndef SPARSESOLV_CORE_SOLVER_CONFIG_HPP
#define SPARSESOLV_CORE_SOLVER_CONFIG_HPP

#include "types.hpp"

namespace sparsesolv {

/**
 * @brief Configuration structure for iterative solvers
 *
 * This struct centralizes all solver parameters that were previously
 * scattered across multiple setter methods. Use this to configure
 * tolerance, iteration limits, preconditioning parameters, and
 * convergence behavior.
 *
 * Example:
 * @code
 * SolverConfig config;
 * config.tolerance = 1e-10;
 * config.max_iterations = 1000;
 * config.shift_parameter = 1.05;  // for IC decomposition
 * @endcode
 */
struct SolverConfig {
    //--------------------------------------------------
    // Convergence criteria
    //--------------------------------------------------

    /// Relative tolerance for convergence (default: 1e-10)
    double tolerance = 1e-10;

    /// Maximum number of iterations (default: 1000)
    int max_iterations = 1000;

    /// Absolute tolerance threshold (solver stops if ||r|| < abs_tolerance)
    double abs_tolerance = 1e-12;

    //--------------------------------------------------
    // Normalization for convergence check
    //--------------------------------------------------

    /// How to normalize the residual for convergence check
    NormType norm_type = NormType::RHS;

    /// Custom normalization value (used when norm_type == Custom)
    double custom_norm = 1.0;

    //--------------------------------------------------
    // Preconditioning parameters
    //--------------------------------------------------

    /// Shift parameter for incomplete Cholesky decomposition
    /// Values > 1.0 improve stability, typical: 1.0-1.2
    /// (Previously called "accera" or acceleration factor)
    double shift_parameter = 1.05;

    /// Enable diagonal scaling (1/sqrt(A[i,i])) before IC factorization
    bool diagonal_scaling = false;

    //--------------------------------------------------
    // Auto-shift parameters for IC decomposition
    //--------------------------------------------------

    /// Enable automatic shift adjustment when IC factorization encounters
    /// small or negative diagonal entries
    bool auto_shift = false;

    /// Increment for automatic shift adjustment (default: 0.01)
    double shift_increment = 0.01;

    /// Maximum allowed shift value during auto-adjustment (default: 5.0)
    double max_shift_value = 5.0;

    /// Threshold below which diagonal is considered too small (triggers shift increase)
    double min_diagonal_threshold = 1e-6;

    /// Replacement value for zero or negative diagonals (default: 1e-10)
    double zero_diagonal_replacement = 1e-10;

    /// Maximum number of shift adjustment trials (default: 100)
    int max_shift_trials = 100;

    //--------------------------------------------------
    // Divergence detection
    //--------------------------------------------------

    /// Strategy for detecting divergence
    DivergenceCheck divergence_check = DivergenceCheck::None;

    /// Multiplier for divergence detection (residual > best * this value triggers count)
    double divergence_threshold = 1000.0;

    /// Number of consecutive iterations above threshold before declaring divergence
    int divergence_count = 100;

    //--------------------------------------------------
    // Result saving options
    //--------------------------------------------------

    /// Save the best result encountered during iteration (useful for non-converging solves)
    bool save_best_result = true;

    /// Save residual history for analysis/debugging
    bool save_residual_history = false;

    //--------------------------------------------------
    // Parallel options
    //--------------------------------------------------

    /// Number of threads for parallel operations (0 = auto)
    int num_threads = 0;

    //--------------------------------------------------
    // ABMC ordering parameters (for parallel IC)
    //--------------------------------------------------

    /// Number of blocks for ABMC ordering
    int abmc_num_blocks = 4;

    /// Number of colors for ABMC ordering
    int abmc_num_colors = 4;

    //--------------------------------------------------
    // Builder pattern for convenient configuration
    //--------------------------------------------------

    SolverConfig& with_tolerance(double tol) {
        tolerance = tol;
        return *this;
    }

    SolverConfig& with_max_iterations(int max_iter) {
        max_iterations = max_iter;
        return *this;
    }

    SolverConfig& with_shift(double shift) {
        shift_parameter = shift;
        return *this;
    }

    SolverConfig& with_diagonal_scaling(bool enable = true) {
        diagonal_scaling = enable;
        return *this;
    }

    SolverConfig& with_residual_history(bool enable = true) {
        save_residual_history = enable;
        return *this;
    }

    SolverConfig& with_divergence_check(DivergenceCheck check, double threshold = 1000.0, int count = 100) {
        divergence_check = check;
        divergence_threshold = threshold;
        divergence_count = count;
        return *this;
    }

    SolverConfig& with_auto_shift(bool enable = true) {
        auto_shift = enable;
        return *this;
    }
};

} // namespace sparsesolv

#endif // SPARSESOLV_CORE_SOLVER_CONFIG_HPP
