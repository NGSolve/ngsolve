/**
 * @file ilu_preconditioner.hpp
 * @brief Incomplete LU (ILU) preconditioner for general matrices
 */

#ifndef SPARSESOLV_PRECONDITIONERS_ILU_PRECONDITIONER_HPP
#define SPARSESOLV_PRECONDITIONERS_ILU_PRECONDITIONER_HPP

#include "../core/preconditioner.hpp"
#include "../core/solver_config.hpp"
#include "../core/level_schedule.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <iostream>
#include <unordered_map>

namespace sparsesolv {

/**
 * @brief Incomplete LU (ILU) preconditioner
 *
 * Computes an incomplete LU factorization: A ≈ L * U
 * where L is a lower triangular matrix with unit diagonal
 * and U is an upper triangular matrix.
 *
 * This preconditioner is suitable for general (non-symmetric) matrices.
 * For symmetric positive definite matrices, IC (Incomplete Cholesky)
 * is more efficient.
 *
 * The shift parameter (α) modifies the diagonal during factorization
 * to improve stability: diag(U) = α * diag(A)
 *
 * Typical values for shift_parameter:
 * - 1.0: Standard ILU(0) - may fail for indefinite matrices
 * - 1.05: Slightly shifted - good default
 * - 1.1-1.2: More stable for difficult problems
 *
 * Example:
 * @code
 * ILUPreconditioner<double> precond(1.05);
 * precond.setup(A);
 *
 * // Apply preconditioner: y = (LU)^{-1} * x
 * precond.apply(x, y, n);
 * @endcode
 *
 * @tparam Scalar The scalar type (double or complex<double>)
 */
template<typename Scalar = double>
class ILUPreconditioner : public Preconditioner<Scalar> {
public:
    /**
     * @brief Construct ILU preconditioner with shift parameter
     * @param shift Shift parameter for diagonal (default: 1.05)
     */
    explicit ILUPreconditioner(double shift = 1.05)
        : shift_parameter_(shift)
    {}

    /**
     * @brief Setup the ILU factorization from matrix A
     *
     * Computes L and U such that A ≈ L * U
     *
     * @param A Sparse matrix (CSR format)
     */
    void setup(const SparseMatrixView<Scalar>& A) override {
        const index_t n = A.rows();
        size_ = n;

        // Copy matrix data - we'll modify it in-place for factorization
        copy_matrix_data(A);

        // Compute ILU factorization
        compute_ilu_factorization();

        // Build level schedules for parallel triangular solves
        // L and U share the same CSR; the schedule builders filter by j < i / j > i
        fwd_schedule_.build_from_lower(row_ptr_.data(), col_idx_.data(), n);
        bwd_schedule_.build_from_upper(row_ptr_.data(), col_idx_.data(), n);

        this->is_setup_ = true;
    }

    /**
     * @brief Apply ILU preconditioner: y = (LU)^{-1} * x
     *
     * Solves L * U * y = x in two steps:
     * 1. Forward substitution: L * z = x
     * 2. Backward substitution: U * y = z
     *
     * @param x Input vector
     * @param y Output vector
     * @param size Vector size
     */
    void apply(const Scalar* x, Scalar* y, index_t size) const override {
        if (!this->is_setup_) {
            throw std::runtime_error("ILUPreconditioner::apply called before setup");
        }

        // Temporary vector for intermediate results
        std::vector<Scalar> z(size);

        // Step 1: Forward substitution - solve L * z = x
        forward_substitution(x, z.data());

        // Step 2: Backward substitution - solve U * y = z
        backward_substitution(z.data(), y);
    }

    std::string name() const override { return "ILU"; }

    /// Get the shift parameter
    double shift_parameter() const { return shift_parameter_; }

    /// Set the shift parameter (must call setup again after changing)
    void set_shift_parameter(double shift) { shift_parameter_ = shift; }

private:
    double shift_parameter_;
    index_t size_ = 0;

    // ILU factorization result stored in modified CSR format
    // The matrix A is stored with L (below diagonal) and U (diagonal and above)
    // L has implicit unit diagonal
    std::vector<index_t> row_ptr_;
    std::vector<index_t> col_idx_;
    std::vector<Scalar> values_;    // Modified to contain L\U factors
    std::vector<index_t> diag_ptr_; // Pointer to diagonal element in each row

    // Level schedules for parallel triangular solves
    LevelSchedule fwd_schedule_;     // For forward substitution (L)
    LevelSchedule bwd_schedule_;     // For backward substitution (U)

    /**
     * @brief Copy matrix data for in-place factorization
     */
    void copy_matrix_data(const SparseMatrixView<Scalar>& A) {
        const index_t n = A.rows();
        const index_t nnz = A.nnz();

        row_ptr_.resize(n + 1);
        col_idx_.resize(nnz);
        values_.resize(nnz);
        diag_ptr_.resize(n);

        // Copy row pointers
        const index_t* src_rp = A.row_ptr();
        parallel_for(n + 1, [&](index_t i) {
            row_ptr_[i] = src_rp[i];
        });

        // Copy column indices and values, find diagonal positions
        const index_t* src_ci = A.col_idx();
        const Scalar* src_v = A.values();
        parallel_for(n, [&](index_t i) {
            index_t diag_pos = row_ptr_[i]; // Default to start of row
            for (index_t k = row_ptr_[i]; k < row_ptr_[i + 1]; ++k) {
                col_idx_[k] = src_ci[k];
                values_[k] = src_v[k];
                if (src_ci[k] == i) {
                    diag_pos = k;
                }
            }
            diag_ptr_[i] = diag_pos;
        });
    }

    /**
     * @brief Compute ILU factorization in-place
     *
     * The factorization is stored in values_ with:
     * - L entries (j < i) in their original positions
     * - U entries (j >= i) in their original positions
     * - L has implicit unit diagonal
     */
    void compute_ilu_factorization() {
        const index_t n = size_;

        // Apply shift to diagonal
        parallel_for(n, [&](index_t i) {
            values_[diag_ptr_[i]] *= static_cast<Scalar>(shift_parameter_);
        });

        // Create a hash map for fast column lookup (per-row)
        // This is sequential due to data dependencies in ILU factorization
        for (index_t i = 0; i < n; ++i) {
            const index_t row_start = row_ptr_[i];
            const index_t row_end = row_ptr_[i + 1];

            // Build column index map for row i
            std::unordered_map<index_t, index_t> col_to_pos;
            for (index_t k = row_start; k < row_end; ++k) {
                col_to_pos[col_idx_[k]] = k;
            }

            // For each column k < i in row i (L part)
            for (index_t kk = row_start; kk < diag_ptr_[i]; ++kk) {
                const index_t k = col_idx_[kk];

                // Get L(i,k) and divide by U(k,k)
                Scalar u_kk = values_[diag_ptr_[k]];
                if (std::abs(u_kk) < std::numeric_limits<double>::epsilon()) {
                    u_kk = std::numeric_limits<double>::epsilon();
                }
                Scalar l_ik = values_[kk] / u_kk;
                values_[kk] = l_ik;

                // Update row i: A(i,j) -= L(i,k) * U(k,j) for j > k
                const index_t k_row_start = row_ptr_[k];
                const index_t k_row_end = row_ptr_[k + 1];

                for (index_t jj = diag_ptr_[k] + 1; jj < k_row_end; ++jj) {
                    const index_t j = col_idx_[jj];
                    auto it = col_to_pos.find(j);
                    if (it != col_to_pos.end()) {
                        // Position exists in row i, update it
                        values_[it->second] -= l_ik * values_[jj];
                    }
                    // If position doesn't exist (fill-in), we drop it (ILU(0))
                }
            }
        }
    }

    /**
     * @brief Forward substitution: solve L * z = x
     * L has unit diagonal (implicit).
     * Uses level scheduling for parallelism.
     */
    void forward_substitution(const Scalar* x, Scalar* z) const {
        for (const auto& level : fwd_schedule_.levels) {
            const index_t level_size = static_cast<index_t>(level.size());
            parallel_for(level_size, [&](index_t idx) {
                const index_t i = level[idx];
                Scalar s = x[i];
                const index_t row_start = row_ptr_[i];
                const index_t diag_pos = diag_ptr_[i];

                // Sum over L(i,j) * z[j] for j < i
                for (index_t k = row_start; k < diag_pos; ++k) {
                    s -= values_[k] * z[col_idx_[k]];
                }

                z[i] = s;  // L has unit diagonal
            });
        }
    }

    /**
     * @brief Backward substitution: solve U * y = z
     * Uses level scheduling for parallelism.
     */
    void backward_substitution(const Scalar* z, Scalar* y) const {
        for (const auto& level : bwd_schedule_.levels) {
            const index_t level_size = static_cast<index_t>(level.size());
            parallel_for(level_size, [&](index_t idx) {
                const index_t i = level[idx];
                Scalar s = z[i];
                const index_t diag_pos = diag_ptr_[i];
                const index_t row_end = row_ptr_[i + 1];

                // Sum over U(i,j) * y[j] for j > i
                for (index_t k = diag_pos + 1; k < row_end; ++k) {
                    s -= values_[k] * y[col_idx_[k]];
                }

                // Divide by diagonal
                y[i] = s / values_[diag_pos];
            });
        }
    }
};

// Type aliases
using ILUPreconditionerD = ILUPreconditioner<double>;
using ILUPreconditionerC = ILUPreconditioner<complex_t>;

} // namespace sparsesolv

#endif // SPARSESOLV_PRECONDITIONERS_ILU_PRECONDITIONER_HPP
