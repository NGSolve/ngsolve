/**
 * @file preconditioner.hpp
 * @brief Base class for preconditioners
 */

#ifndef SPARSESOLV_CORE_PRECONDITIONER_HPP
#define SPARSESOLV_CORE_PRECONDITIONER_HPP

#include "types.hpp"
#include "constants.hpp"
#include "parallel.hpp"
#include "sparse_matrix_view.hpp"
#include <string>
#include <memory>

namespace sparsesolv {

/**
 * @brief Abstract base class for preconditioners
 *
 * A preconditioner approximates M^{-1} for a matrix A, where M is chosen
 * such that M^{-1}A has better spectral properties than A. This accelerates
 * convergence of iterative solvers.
 *
 * To implement a preconditioner:
 * 1. Override setup() to compute the preconditioner from matrix A
 * 2. Override apply() to compute y = M^{-1} * x
 * 3. Override name() to return a descriptive name
 *
 * Example:
 * @code
 * class MyPreconditioner : public Preconditioner<double> {
 * public:
 *     void setup(const SparseMatrixView<double>& A) override {
 *         // Compute preconditioner from A
 *     }
 *
 *     void apply(const double* x, double* y, index_t size) const override {
 *         // Compute y = M^{-1} * x
 *     }
 *
 *     std::string name() const override { return "MyPreconditioner"; }
 * };
 * @endcode
 *
 * @tparam Scalar The scalar type (double or complex<double>)
 */
template<typename Scalar = double>
class Preconditioner {
public:
    using value_type = Scalar;

    virtual ~Preconditioner() = default;

    /**
     * @brief Setup the preconditioner from matrix A
     *
     * This method computes any factorization or decomposition needed
     * for the preconditioner. It must be called before apply().
     *
     * @param A The matrix to precondition
     */
    virtual void setup(const SparseMatrixView<Scalar>& A) = 0;

    /**
     * @brief Apply the preconditioner: y = M^{-1} * x
     *
     * @param x Input vector
     * @param y Output vector
     * @param size Vector size
     */
    virtual void apply(const Scalar* x, Scalar* y, index_t size) const = 0;

    /**
     * @brief Apply the preconditioner with std::vector
     */
    void apply(const std::vector<Scalar>& x, std::vector<Scalar>& y) const {
        assert(x.size() == y.size());
        apply(x.data(), y.data(), static_cast<index_t>(x.size()));
    }

    /**
     * @brief Get the name of the preconditioner
     */
    virtual std::string name() const = 0;

    /**
     * @brief Check if the preconditioner has been set up
     */
    virtual bool is_ready() const { return is_setup_; }

protected:
    bool is_setup_ = false;
};

/**
 * @brief Identity preconditioner (no preconditioning)
 *
 * This is a null preconditioner that simply copies the input to output.
 * Useful as a baseline for comparison.
 */
template<typename Scalar = double>
class IdentityPreconditioner : public Preconditioner<Scalar> {
public:
    void setup(const SparseMatrixView<Scalar>& /*A*/) override {
        this->is_setup_ = true;
    }

    void apply(const Scalar* x, Scalar* y, index_t size) const override {
        std::copy(x, x + size, y);
    }

    std::string name() const override { return "Identity"; }
};

/**
 * @brief Jacobi (diagonal) preconditioner
 *
 * M = diag(A), so M^{-1} * x = x ./ diag(A)
 */
template<typename Scalar = double>
class JacobiPreconditioner : public Preconditioner<Scalar> {
public:
    void setup(const SparseMatrixView<Scalar>& A) override {
        index_t n = A.rows();
        inv_diag_.resize(n);
        for (index_t i = 0; i < n; ++i) {
            Scalar d = A.diagonal(i);
            // Avoid division by zero
            inv_diag_[i] = (std::abs(d) > constants::MIN_DIAGONAL_TOLERANCE) ? Scalar(1) / d : Scalar(1);
        }
        this->is_setup_ = true;
    }

    void apply(const Scalar* x, Scalar* y, index_t size) const override {
        parallel_for(size, [&](index_t i) {
            y[i] = inv_diag_[i] * x[i];
        });
    }

    std::string name() const override { return "Jacobi"; }

private:
    std::vector<Scalar> inv_diag_;
};

// Type aliases
using PreconditionerD = Preconditioner<double>;
using PreconditionerC = Preconditioner<complex_t>;

} // namespace sparsesolv

#endif // SPARSESOLV_CORE_PRECONDITIONER_HPP
