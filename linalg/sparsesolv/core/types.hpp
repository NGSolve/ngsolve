/**
 * @file types.hpp
 * @brief Basic type definitions for SparseSolv library
 */

#ifndef SPARSESOLV_CORE_TYPES_HPP
#define SPARSESOLV_CORE_TYPES_HPP

#include <cstdint>
#include <complex>
#include <vector>
#include <string>

namespace sparsesolv {

// Index type for sparse matrix indices
using index_t = std::int32_t;

// Complex type alias
using complex_t = std::complex<double>;

/**
 * @brief Result of an iterative solver
 */
struct SolverResult {
    bool converged = false;           ///< Whether the solver converged
    int iterations = 0;               ///< Number of iterations performed
    double final_residual = 0.0;      ///< Final relative residual
    std::vector<double> residual_history;  ///< Residual at each iteration (optional)

    /// Check if the solve was successful
    explicit operator bool() const { return converged; }
};

/**
 * @brief Type of norm used for convergence check
 */
enum class NormType {
    RHS,              ///< Normalize by ||b|| (right-hand side norm)
    InitialResidual,  ///< Normalize by ||r_0|| (initial residual)
    Custom            ///< Use a custom normalization value
};

/**
 * @brief Divergence detection strategy
 */
enum class DivergenceCheck {
    None,             ///< No divergence check, run until max iterations
    StagnationCount   ///< Stop if residual stagnates for too many iterations
};

} // namespace sparsesolv

#endif // SPARSESOLV_CORE_TYPES_HPP
