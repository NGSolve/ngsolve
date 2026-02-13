/**
 * @file level_schedule.hpp
 * @brief Level scheduling for parallel triangular solves
 *
 * Computes dependency levels for rows of a triangular matrix so that
 * rows at the same level can be processed in parallel. This enables
 * parallel forward/backward substitution using parallel_for.
 *
 * Algorithm:
 *   For lower triangular L:
 *     level[i] = max(level[j] for j in L[i,:] where j < i) + 1
 *   Rows at the same level are independent and can run in parallel.
 *
 * For upper triangular (backward solve), build levels from the transpose
 * pattern, processing from the last row backward.
 */

#ifndef SPARSESOLV_CORE_LEVEL_SCHEDULE_HPP
#define SPARSESOLV_CORE_LEVEL_SCHEDULE_HPP

#include "types.hpp"
#include "parallel.hpp"
#include <vector>
#include <algorithm>

namespace sparsesolv {

/**
 * @brief Level schedule for parallel triangular solves
 *
 * Stores rows grouped by dependency level. Rows within the same level
 * are independent and can be processed in parallel via parallel_for.
 *
 * Construction cost: O(nnz), done once per factorization.
 * Usage: iterate levels sequentially, within each level use parallel_for.
 */
struct LevelSchedule {
    /// levels[k] = vector of row indices that can be processed in parallel at step k
    std::vector<std::vector<index_t>> levels;

    /// Total number of rows
    index_t size = 0;

    /**
     * @brief Build level schedule from lower triangular CSR pattern
     *
     * For forward substitution: L * y = b, processing rows 0..n-1.
     * Row i depends on rows j < i where L[i,j] != 0.
     *
     * @param row_ptr CSR row pointer array (size n+1)
     * @param col_idx CSR column index array
     * @param n Number of rows
     */
    void build_from_lower(const index_t* row_ptr, const index_t* col_idx, index_t n) {
        size = n;
        if (n <= 0) {
            levels.clear();
            return;
        }

        // Compute level for each row
        std::vector<index_t> row_level(n, 0);
        index_t max_level = 0;

        for (index_t i = 0; i < n; ++i) {
            index_t max_dep = -1;
            // Scan off-diagonal entries (j < i) in row i
            for (index_t k = row_ptr[i]; k < row_ptr[i + 1]; ++k) {
                index_t j = col_idx[k];
                if (j < i) {
                    max_dep = std::max(max_dep, row_level[j]);
                }
            }
            row_level[i] = max_dep + 1;
            max_level = std::max(max_level, row_level[i]);
        }

        // Group rows by level
        levels.resize(max_level + 1);
        for (auto& lev : levels) lev.clear();

        for (index_t i = 0; i < n; ++i) {
            levels[row_level[i]].push_back(i);
        }
    }

    /**
     * @brief Build level schedule from upper triangular CSR pattern
     *
     * For backward substitution: U * y = b, processing rows n-1..0.
     * Row i depends on rows j > i where U[i,j] != 0.
     *
     * The levels are ordered so that level 0 should be processed first
     * (contains the last rows), then level 1, etc.
     *
     * @param row_ptr CSR row pointer array (size n+1)
     * @param col_idx CSR column index array
     * @param n Number of rows
     */
    void build_from_upper(const index_t* row_ptr, const index_t* col_idx, index_t n) {
        size = n;
        if (n <= 0) {
            levels.clear();
            return;
        }

        // Compute level for each row (processing from last to first)
        std::vector<index_t> row_level(n, 0);
        index_t max_level = 0;

        for (index_t i = n; i-- > 0;) {
            index_t max_dep = -1;
            // Scan off-diagonal entries (j > i) in row i
            for (index_t k = row_ptr[i]; k < row_ptr[i + 1]; ++k) {
                index_t j = col_idx[k];
                if (j > i) {
                    max_dep = std::max(max_dep, row_level[j]);
                }
            }
            row_level[i] = max_dep + 1;
            max_level = std::max(max_level, row_level[i]);
        }

        // Group rows by level
        levels.resize(max_level + 1);
        for (auto& lev : levels) lev.clear();

        for (index_t i = 0; i < n; ++i) {
            levels[row_level[i]].push_back(i);
        }
    }

    /// Check if the schedule has been built
    bool is_built() const { return !levels.empty(); }

    /// Number of levels (= sequential steps needed)
    index_t num_levels() const { return static_cast<index_t>(levels.size()); }
};

} // namespace sparsesolv

#endif // SPARSESOLV_CORE_LEVEL_SCHEDULE_HPP
