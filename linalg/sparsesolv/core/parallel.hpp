/**
 * @file parallel.hpp
 * @brief Portable parallel primitives for SparseSolv
 *
 * Provides compile-time dispatch to:
 * - NGSolve TaskManager (when SPARSESOLV_USE_NGSOLVE_TASKMANAGER is defined)
 * - OpenMP (when _OPENMP is defined)
 * - Serial fallback (otherwise)
 *
 * This allows SparseSolv to work both as a standalone library and
 * integrated into NGSolve without OpenMP/TaskManager thread pool conflicts.
 */

#ifndef SPARSESOLV_CORE_PARALLEL_HPP
#define SPARSESOLV_CORE_PARALLEL_HPP

#include "types.hpp"
#include <type_traits>
#include <complex>

#ifdef SPARSESOLV_USE_NGSOLVE_TASKMANAGER
  #include <core/taskmanager.hpp>
#elif defined(_OPENMP)
  #include <omp.h>
#endif

namespace sparsesolv {

/**
 * @brief Parallel for loop over [0, n)
 *
 * Dispatches to ngcore::ParallelFor, OpenMP, or serial loop.
 */
template<typename FUNC>
inline void parallel_for(index_t n, FUNC f) {
#ifdef SPARSESOLV_USE_NGSOLVE_TASKMANAGER
    if (n > 0)
        ngcore::ParallelFor(static_cast<size_t>(n),
            [&](size_t i) { f(static_cast<index_t>(i)); });
#elif defined(_OPENMP)
    #pragma omp parallel for schedule(static)
    for (index_t i = 0; i < n; ++i)
        f(i);
#else
    for (index_t i = 0; i < n; ++i)
        f(i);
#endif
}

/**
 * @brief Parallel reduction (summation) over [0, n)
 *
 * Computes sum of f(i) for i in [0, n).
 * Uses ngcore::ParallelReduce (supports all types including complex),
 * OpenMP reduction (double only on MSVC), or serial loop.
 *
 * @tparam T Result type (double, complex<double>, etc.)
 * @tparam FUNC Function type: index_t -> T
 * @param n Loop bound
 * @param f Function to evaluate at each index
 * @param init Initial value for reduction
 * @return Sum of f(i) for i in [0, n)
 */
template<typename T, typename FUNC>
inline T parallel_reduce_sum(index_t n, FUNC f, T init = T(0)) {
#ifdef SPARSESOLV_USE_NGSOLVE_TASKMANAGER
    if (n <= 0) return init;
    return ngcore::ParallelReduce(static_cast<size_t>(n),
        [&](size_t i) { return f(static_cast<index_t>(i)); },
        [](T a, T b) { return a + b; }, init);
#elif defined(_OPENMP)
    T sum = init;
    // MSVC OpenMP (v2.0) only supports reduction on scalar arithmetic types
    if constexpr (std::is_same_v<T, double>) {
        #pragma omp parallel for reduction(+:sum)
        for (index_t i = 0; i < n; ++i)
            sum += f(i);
    } else {
        // Complex or other types: serial fallback on MSVC
        for (index_t i = 0; i < n; ++i)
            sum += f(i);
    }
    return sum;
#else
    T sum = init;
    for (index_t i = 0; i < n; ++i)
        sum += f(i);
    return sum;
#endif
}

} // namespace sparsesolv

#endif // SPARSESOLV_CORE_PARALLEL_HPP
