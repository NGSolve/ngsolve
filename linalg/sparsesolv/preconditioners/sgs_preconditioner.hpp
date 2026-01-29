/**
 * @file sgs_preconditioner.hpp
 * @brief Symmetric Gauss-Seidel (SGS) preconditioner
 */

#ifndef SPARSESOLV_PRECONDITIONERS_SGS_PRECONDITIONER_HPP
#define SPARSESOLV_PRECONDITIONERS_SGS_PRECONDITIONER_HPP

#include "../core/preconditioner.hpp"
#include <algorithm>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace sparsesolv {

/**
 * @brief Symmetric Gauss-Seidel (SGS) preconditioner
 *
 * Computes the SGS preconditioner: M = (D + L) * D^{-1} * (D + U)
 * where A = L + D + U (lower, diagonal, upper decomposition).
 *
 * The application M^{-1} * x is computed as:
 * 1. Forward sweep: solve (D + L) * y = x
 * 2. Diagonal scaling: z = D * y
 * 3. Backward sweep: solve (D + U) * result = z
 *
 * SGS is a good general-purpose preconditioner that:
 * - Requires no additional fill-in
 * - Works well with MRTR (SGS-MRTR)
 * - Is symmetric, so can be used with CG
 *
 * @note This class owns a copy of the matrix data, so the original matrix
 *       can be safely destroyed after setup().
 * @note This class is NOT thread-safe for concurrent apply() calls.
 *       Use separate preconditioner instances for each thread.
 *
 * @tparam Scalar The scalar type (double or complex<double>)
 */
template<typename Scalar = double>
class SGSPreconditioner : public Preconditioner<Scalar> {
public:
    SGSPreconditioner() = default;

    // Disable copy (owns resources)
    SGSPreconditioner(const SGSPreconditioner&) = delete;
    SGSPreconditioner& operator=(const SGSPreconditioner&) = delete;

    // Enable move
    SGSPreconditioner(SGSPreconditioner&&) noexcept = default;
    SGSPreconditioner& operator=(SGSPreconditioner&&) noexcept = default;

    /**
     * @brief Setup the SGS preconditioner from matrix A
     *
     * Copies the matrix data internally for safe storage.
     *
     * @param A Sparse matrix (CSR format)
     */
    void setup(const SparseMatrixView<Scalar>& A) override {
        const index_t n = A.rows();
        size_ = n;

        // Copy CSR data (owns the data to avoid dangling pointers)
        const index_t nnz = A.nnz();
        row_ptr_.resize(n + 1);
        col_idx_.resize(nnz);
        values_.resize(nnz);

        // Copy row pointers (OpenMP parallel)
        const index_t* src_rp = A.row_ptr();
        #pragma omp parallel for schedule(static)
        for (index_t i = 0; i <= n; ++i) {
            row_ptr_[i] = src_rp[i];
        }

        // Copy column indices and values (OpenMP parallel)
        const index_t* src_ci = A.col_idx();
        const Scalar* src_v = A.values();
        #pragma omp parallel for schedule(static)
        for (index_t i = 0; i < nnz; ++i) {
            col_idx_[i] = src_ci[i];
            values_[i] = src_v[i];
        }

        // Extract and invert diagonal (OpenMP parallel)
        inv_diag_.resize(n);
        #pragma omp parallel for schedule(static)
        for (index_t i = 0; i < n; ++i) {
            Scalar d = A.diagonal(i);
            if (std::abs(d) > 1e-15) {
                inv_diag_[i] = Scalar(1) / d;
            } else {
                inv_diag_[i] = Scalar(1);
            }
        }

        this->is_setup_ = true;
    }

    /**
     * @brief Apply SGS preconditioner: y = M^{-1} * x
     *
     * M^{-1} = (D + U)^{-1} * D * (D + L)^{-1}
     *
     * @param x Input vector
     * @param y Output vector
     * @param size Vector size
     *
     * @note NOT thread-safe. Use separate instances for concurrent access.
     */
    void apply(const Scalar* x, Scalar* y, index_t size) const override {
        if (!this->is_setup_) {
            throw std::runtime_error("SGSPreconditioner::apply called before setup");
        }

        // Allocate temporary vectors for this call (thread-safe)
        std::vector<Scalar> temp(size_);
        std::vector<Scalar> scaled(size_);

        // Step 1: Forward Gauss-Seidel sweep (sequential - data dependency)
        // Solve (D + L) * temp = x  =>  temp = (D + L)^{-1} * x
        forward_sweep(x, temp.data());

        // Step 2: Diagonal scaling (OpenMP parallel - no data dependency)
        // scaled = D * temp
        #pragma omp parallel for schedule(static)
        for (index_t i = 0; i < size_; ++i) {
            // inv_diag_[i] = 1/D[i,i], so D[i,i] = 1/inv_diag_[i]
            scaled[i] = temp[i] / inv_diag_[i];
        }

        // Step 3: Backward Gauss-Seidel sweep (sequential - data dependency)
        // Solve (D + U) * y = scaled  =>  y = (D + U)^{-1} * D * (D + L)^{-1} * x
        backward_sweep(scaled.data(), y);
    }

    std::string name() const override { return "SGS"; }

private:
    index_t size_ = 0;

    // Owned CSR data (copied from input matrix)
    std::vector<index_t> row_ptr_;
    std::vector<index_t> col_idx_;
    std::vector<Scalar> values_;

    std::vector<Scalar> inv_diag_;

    /**
     * @brief Forward Gauss-Seidel sweep: solve (D + L) * y = x
     *
     * For i = 0, 1, ..., n-1:
     *   y[i] = (x[i] - sum_{j<i} A[i,j] * y[j]) / A[i,i]
     */
    void forward_sweep(const Scalar* x, Scalar* y) const {
        const index_t* rp = row_ptr_.data();
        const index_t* ci = col_idx_.data();
        const Scalar* v = values_.data();

        for (index_t i = 0; i < size_; ++i) {
            Scalar sum = x[i];
            for (index_t k = rp[i]; k < rp[i + 1]; ++k) {
                index_t j = ci[k];
                if (j < i) {
                    sum -= v[k] * y[j];
                }
            }
            y[i] = sum * inv_diag_[i];
        }
    }

    /**
     * @brief Backward Gauss-Seidel sweep: solve (D + U) * y = x
     *
     * For i = n-1, n-2, ..., 0:
     *   y[i] = (x[i] - sum_{j>i} A[i,j] * y[j]) / A[i,i]
     */
    void backward_sweep(const Scalar* x, Scalar* y) const {
        const index_t* rp = row_ptr_.data();
        const index_t* ci = col_idx_.data();
        const Scalar* v = values_.data();

        for (index_t i = size_; i-- > 0;) {
            Scalar sum = x[i];
            for (index_t k = rp[i]; k < rp[i + 1]; ++k) {
                index_t j = ci[k];
                if (j > i) {
                    sum -= v[k] * y[j];
                }
            }
            y[i] = sum * inv_diag_[i];
        }
    }
};

// Type aliases
using SGSPreconditionerD = SGSPreconditioner<double>;
using SGSPreconditionerC = SGSPreconditioner<complex_t>;

} // namespace sparsesolv

#endif // SPARSESOLV_PRECONDITIONERS_SGS_PRECONDITIONER_HPP
