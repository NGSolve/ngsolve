/**
 * @file sparse_matrix_view.hpp
 * @brief Zero-copy view into sparse matrices (CSR format)
 */

#ifndef SPARSESOLV_CORE_SPARSE_MATRIX_VIEW_HPP
#define SPARSESOLV_CORE_SPARSE_MATRIX_VIEW_HPP

#include "types.hpp"
#include "parallel.hpp"
#include <cassert>
#include <stdexcept>

namespace sparsesolv {

/**
 * @brief Non-owning view into a sparse matrix in CSR format
 *
 * This class provides a unified interface to sparse matrices without
 * copying data. It can wrap matrices from various sources:
 * - Raw CSR arrays (row_ptr, col_idx, values)
 * - scipy.sparse.csr_matrix (via pybind11)
 * - Eigen::SparseMatrix
 * - NGSolve's SparseMatrix
 *
 * The CSR format stores:
 * - row_ptr[i]: Index into col_idx/values where row i starts
 * - col_idx[j]: Column index for entry j
 * - values[j]: Value for entry j
 *
 * Example:
 * @code
 * // Wrap existing CSR data
 * SparseMatrixView<double> view(n, n, row_ptr, col_idx, values);
 *
 * // Matrix-vector product: y = A * x
 * view.multiply(x, y);
 *
 * // Access elements
 * double aij = view(i, j);  // O(log nnz_per_row) lookup
 * @endcode
 *
 * @warning LIFETIME MANAGEMENT: This is a non-owning view. The underlying
 *          arrays (row_ptr, col_idx, values) MUST remain valid for the
 *          entire lifetime of the view. If the source data is deallocated
 *          or reallocated while the view exists, the behavior is undefined.
 *
 * @note For preconditioners that need to outlive the source matrix,
 *       use ICPreconditioner or SGSPreconditioner, which copy the data.
 *
 * @tparam Scalar The scalar type (double or complex<double>)
 */
template<typename Scalar = double>
class SparseMatrixView {
public:
    using value_type = Scalar;
    using index_type = index_t;

    /**
     * @brief Construct a view from CSR arrays
     * @param rows Number of rows
     * @param cols Number of columns
     * @param row_ptr Row pointer array (size: rows + 1)
     * @param col_idx Column index array (size: nnz)
     * @param values Value array (size: nnz)
     */
    SparseMatrixView(index_t rows, index_t cols,
                     const index_t* row_ptr,
                     const index_t* col_idx,
                     const Scalar* values)
        : rows_(rows)
        , cols_(cols)
        , row_ptr_(row_ptr)
        , col_idx_(col_idx)
        , values_(values)
    {
        assert(row_ptr != nullptr);
        assert(col_idx != nullptr);
        assert(values != nullptr);
    }

    /// Default constructor (creates an empty view)
    SparseMatrixView()
        : rows_(0), cols_(0), row_ptr_(nullptr), col_idx_(nullptr), values_(nullptr)
    {}

    // Accessors

    /// Number of rows
    index_t rows() const { return rows_; }

    /// Number of columns
    index_t cols() const { return cols_; }

    /// Number of non-zeros
    index_t nnz() const { return row_ptr_ ? row_ptr_[rows_] : 0; }

    /// Check if the view is valid
    bool is_valid() const { return row_ptr_ != nullptr; }

    /// Raw access to row pointers
    const index_t* row_ptr() const { return row_ptr_; }

    /// Raw access to column indices
    const index_t* col_idx() const { return col_idx_; }

    /// Raw access to values
    const Scalar* values() const { return values_; }

    /**
     * @brief Get the start and end indices for row i
     * @param i Row index
     * @return Pair of (start, end) indices into col_idx/values
     */
    std::pair<index_t, index_t> row_range(index_t i) const {
        assert(i >= 0 && i < rows_);
        return {row_ptr_[i], row_ptr_[i + 1]};
    }

    /**
     * @brief Access element (i, j) - O(log nnz_per_row) binary search
     * @param i Row index
     * @param j Column index
     * @return Value at (i, j), or 0 if not present
     */
    Scalar operator()(index_t i, index_t j) const {
        auto [start, end] = row_range(i);
        // Binary search for column j
        auto it = std::lower_bound(col_idx_ + start, col_idx_ + end, j);
        if (it != col_idx_ + end && *it == j) {
            return values_[it - col_idx_];
        }
        return Scalar(0);
    }

    /**
     * @brief Matrix-vector product: y = A * x
     * @param x Input vector (size: cols)
     * @param y Output vector (size: rows)
     */
    void multiply(const Scalar* x, Scalar* y) const {
        parallel_for(rows_, [&](index_t i) {
            Scalar sum = Scalar(0);
            for (index_t k = row_ptr_[i]; k < row_ptr_[i + 1]; ++k) {
                sum += values_[k] * x[col_idx_[k]];
            }
            y[i] = sum;
        });
    }

    /**
     * @brief Matrix-vector product: y = A * x (with std::vector)
     */
    void multiply(const std::vector<Scalar>& x, std::vector<Scalar>& y) const {
        assert(static_cast<index_t>(x.size()) >= cols_);
        assert(static_cast<index_t>(y.size()) >= rows_);
        multiply(x.data(), y.data());
    }

    /**
     * @brief Get diagonal element
     * @param i Row/column index
     * @return Diagonal value A(i,i)
     */
    Scalar diagonal(index_t i) const {
        return (*this)(i, i);
    }

    /**
     * @brief Extract diagonal as a vector
     * @param diag Output vector (size: min(rows, cols))
     */
    void get_diagonal(Scalar* diag) const {
        index_t n = std::min(rows_, cols_);
        for (index_t i = 0; i < n; ++i) {
            diag[i] = diagonal(i);
        }
    }

private:
    index_t rows_;
    index_t cols_;
    const index_t* row_ptr_;
    const index_t* col_idx_;
    const Scalar* values_;
};

// Type aliases
using SparseMatrixViewD = SparseMatrixView<double>;
using SparseMatrixViewC = SparseMatrixView<complex_t>;

} // namespace sparsesolv

#endif // SPARSESOLV_CORE_SPARSE_MATRIX_VIEW_HPP
