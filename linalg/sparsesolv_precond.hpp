/**
 * @file sparsesolv_precond.hpp
 * @brief SparseSolv preconditioners integrated with NGSolve
 *
 * Provides IC (Incomplete Cholesky), ILU (Incomplete LU), and
 * SGS (Symmetric Gauss-Seidel) preconditioners for use with NGSolve's Krylov solvers.
 *
 * These preconditioners support Dirichlet boundary conditions through
 * the freedofs parameter (BitArray).
 *
 * @code
 * from ngsolve import *
 * from ngsolve.krylovspace import CGSolver
 *
 * # Create bilinear form with Dirichlet BC
 * fes = H1(mesh, order=2, dirichlet="left|right|top|bottom")
 * a = BilinearForm(fes)
 * a += grad(u)*grad(v)*dx
 * a.Assemble()
 *
 * # Create IC preconditioner with FreeDofs
 * pre = ICPreconditioner(a.mat, freedofs=fes.FreeDofs(), shift=1.05)
 * pre.Update()
 *
 * # Use with CGSolver
 * inv = CGSolver(a.mat, pre, printrates=True, tol=1e-10)
 * gfu.vec.data = inv * f.vec
 * @endcode
 *
 * Based on JP-MARs/SparseSolv (https://github.com/JP-MARs/SparseSolv)
 */

#ifndef NGSOLVE_SPARSESOLV_PRECOND_HPP
#define NGSOLVE_SPARSESOLV_PRECOND_HPP

#include <basematrix.hpp>
#include <sparsematrix.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

// SparseSolv headers
#include "sparsesolv/sparsesolv.hpp"

namespace ngla {

// ============================================================================
// Helper Functions
// ============================================================================

/**
 * @brief Helper to extract raw data pointer from NGSolve vectors
 */
template<typename SCAL>
const SCAL* GetVectorData(const BaseVector& vec) {
    return vec.FV<SCAL>().Data();
}

template<typename SCAL>
SCAL* GetVectorData(BaseVector& vec) {
    return vec.FV<SCAL>().Data();
}

// ============================================================================
// IC Preconditioner with FreeDofs support
// ============================================================================

/**
 * @brief Incomplete Cholesky Preconditioner for NGSolve
 *
 * Uses shifted IC decomposition for improved stability.
 * Supports Dirichlet boundary conditions through freedofs parameter.
 *
 * @tparam SCAL Scalar type (double or Complex)
 */
template<typename SCAL = double>
class SparseSolvICPreconditioner : public BaseMatrix {
public:
    /**
     * @brief Construct IC preconditioner from NGSolve sparse matrix
     * @param mat The sparse matrix to precondition
     * @param freedofs BitArray indicating free DOFs (nullptr = all DOFs are free)
     * @param shift Shift parameter for IC decomposition (default: 1.05)
     */
    SparseSolvICPreconditioner(shared_ptr<SparseMatrix<SCAL>> mat,
                                shared_ptr<BitArray> freedofs = nullptr,
                                double shift = 1.05)
        : mat_(mat)
        , freedofs_(freedofs ? make_shared<BitArray>(*freedofs) : nullptr)
        , shift_(shift)
        , height_(static_cast<sparsesolv::index_t>(mat->Height()))
        , width_(static_cast<sparsesolv::index_t>(mat->Width()))
        , precond_(std::make_shared<sparsesolv::ICPreconditioner<SCAL>>(shift))
    {}

    /**
     * @brief Update the preconditioner (recompute factorization)
     */
    void Update() {
        const auto& firsti = mat_->GetFirstArray();
        const auto& colnr = mat_->GetColIndices();
        const auto& values = mat_->GetValues();

        const size_t nrows = firsti.Size();
        const size_t nnz = colnr.Size();

        // Convert row pointers (OpenMP parallel)
        row_ptr_.resize(nrows);
        #pragma omp parallel for schedule(static)
        for (ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(nrows); ++i) {
            row_ptr_[i] = static_cast<sparsesolv::index_t>(firsti[i]);
        }

        // Convert column indices (OpenMP parallel)
        col_idx_.resize(nnz);
        #pragma omp parallel for schedule(static)
        for (ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(nnz); ++i) {
            col_idx_[i] = static_cast<sparsesolv::index_t>(colnr[i]);
        }

        // Create modified values array if freedofs is specified
        // For constrained DOFs, set diagonal to 1 and off-diagonals to 0
        // Also zero out column entries pointing to constrained DOFs (for symmetry)
        if (freedofs_) {
            modified_values_.resize(values.Size());

            // Copy values (OpenMP parallel)
            #pragma omp parallel for schedule(static)
            for (ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(values.Size()); ++i) {
                modified_values_[i] = values[i];
            }

            // Modify matrix entries for constrained DOFs (OpenMP parallel by row)
            #pragma omp parallel for schedule(dynamic, 64)
            for (sparsesolv::index_t i = 0; i < height_; ++i) {
                for (sparsesolv::index_t k = row_ptr_[i]; k < row_ptr_[i + 1]; ++k) {
                    sparsesolv::index_t j = col_idx_[k];
                    if (!freedofs_->Test(i)) {
                        // Constrained row: set to identity row
                        if (j == i) {
                            modified_values_[k] = SCAL(1);  // Diagonal = 1
                        } else {
                            modified_values_[k] = SCAL(0);  // Off-diagonal = 0
                        }
                    } else if (!freedofs_->Test(j)) {
                        // Free row but column points to constrained DOF: zero it
                        modified_values_[k] = SCAL(0);
                    }
                }
            }

            sparsesolv::SparseMatrixView<SCAL> view(
                height_, width_,
                row_ptr_.data(),
                col_idx_.data(),
                modified_values_.data()
            );
            precond_->setup(view);
        } else {
            sparsesolv::SparseMatrixView<SCAL> view(
                height_, width_,
                row_ptr_.data(),
                col_idx_.data(),
                values.Data()
            );
            precond_->setup(view);
        }
    }

    /**
     * @brief Apply preconditioner: y = M^{-1} * x
     */
    void Mult(const BaseVector& x, BaseVector& y) const override {
        const SCAL* x_data = GetVectorData<SCAL>(x);
        SCAL* y_data = GetVectorData<SCAL>(y);

        if (freedofs_) {
            // Apply preconditioner
            precond_->apply(x_data, y_data, height_);

            // Zero out constrained DOFs (OpenMP parallel)
            #pragma omp parallel for schedule(static)
            for (sparsesolv::index_t i = 0; i < height_; ++i) {
                if (!freedofs_->Test(i)) {
                    y_data[i] = SCAL(0);
                }
            }
        } else {
            precond_->apply(x_data, y_data, height_);
        }
    }

    void MultAdd(double s, const BaseVector& x, BaseVector& y) const override {
        const SCAL* x_data = GetVectorData<SCAL>(x);
        SCAL* y_data = GetVectorData<SCAL>(y);

        std::vector<SCAL> temp(height_);
        precond_->apply(x_data, temp.data(), height_);

        if (freedofs_) {
            #pragma omp parallel for schedule(static)
            for (sparsesolv::index_t i = 0; i < height_; ++i) {
                if (freedofs_->Test(i)) {
                    y_data[i] += s * temp[i];
                }
            }
        } else {
            #pragma omp parallel for schedule(static)
            for (sparsesolv::index_t i = 0; i < height_; ++i) {
                y_data[i] += s * temp[i];
            }
        }
    }

    int VHeight() const override { return static_cast<int>(height_); }
    int VWidth() const override { return static_cast<int>(width_); }

    AutoVector CreateRowVector() const override {
        return mat_->CreateRowVector();
    }

    AutoVector CreateColVector() const override {
        return mat_->CreateColVector();
    }

    double GetShift() const { return shift_; }

    void SetShift(double shift) {
        shift_ = shift;
        precond_->set_shift_parameter(shift);
    }

private:
    shared_ptr<SparseMatrix<SCAL>> mat_;
    shared_ptr<BitArray> freedofs_;
    double shift_;
    sparsesolv::index_t height_;
    sparsesolv::index_t width_;
    shared_ptr<sparsesolv::ICPreconditioner<SCAL>> precond_;

    std::vector<sparsesolv::index_t> row_ptr_;
    std::vector<sparsesolv::index_t> col_idx_;
    std::vector<SCAL> modified_values_;
};

// ============================================================================
// SGS Preconditioner with FreeDofs support
// ============================================================================

/**
 * @brief Symmetric Gauss-Seidel Preconditioner for NGSolve
 *
 * Supports Dirichlet boundary conditions through freedofs parameter.
 *
 * @tparam SCAL Scalar type (double or Complex)
 */
template<typename SCAL = double>
class SparseSolvSGSPreconditioner : public BaseMatrix {
public:
    /**
     * @brief Construct SGS preconditioner from NGSolve sparse matrix
     * @param mat The sparse matrix to precondition
     * @param freedofs BitArray indicating free DOFs (nullptr = all DOFs are free)
     */
    SparseSolvSGSPreconditioner(shared_ptr<SparseMatrix<SCAL>> mat,
                                 shared_ptr<BitArray> freedofs = nullptr)
        : mat_(mat)
        , freedofs_(freedofs ? make_shared<BitArray>(*freedofs) : nullptr)
        , height_(static_cast<sparsesolv::index_t>(mat->Height()))
        , width_(static_cast<sparsesolv::index_t>(mat->Width()))
        , precond_(std::make_shared<sparsesolv::SGSPreconditioner<SCAL>>())
    {}

    /**
     * @brief Update the preconditioner
     */
    void Update() {
        const auto& firsti = mat_->GetFirstArray();
        const auto& colnr = mat_->GetColIndices();
        const auto& values = mat_->GetValues();

        const size_t nrows = firsti.Size();
        const size_t nnz = colnr.Size();

        // Convert row pointers (OpenMP parallel)
        row_ptr_.resize(nrows);
        #pragma omp parallel for schedule(static)
        for (ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(nrows); ++i) {
            row_ptr_[i] = static_cast<sparsesolv::index_t>(firsti[i]);
        }

        // Convert column indices (OpenMP parallel)
        col_idx_.resize(nnz);
        #pragma omp parallel for schedule(static)
        for (ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(nnz); ++i) {
            col_idx_[i] = static_cast<sparsesolv::index_t>(colnr[i]);
        }

        // Create modified values array if freedofs is specified
        // For constrained DOFs, set diagonal to 1 and off-diagonals to 0
        // Also zero out column entries pointing to constrained DOFs (for symmetry)
        if (freedofs_) {
            modified_values_.resize(values.Size());

            // Copy values (OpenMP parallel)
            #pragma omp parallel for schedule(static)
            for (ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(values.Size()); ++i) {
                modified_values_[i] = values[i];
            }

            // Modify matrix entries for constrained DOFs (OpenMP parallel by row)
            #pragma omp parallel for schedule(dynamic, 64)
            for (sparsesolv::index_t i = 0; i < height_; ++i) {
                for (sparsesolv::index_t k = row_ptr_[i]; k < row_ptr_[i + 1]; ++k) {
                    sparsesolv::index_t j = col_idx_[k];
                    if (!freedofs_->Test(i)) {
                        // Constrained row: set to identity row
                        if (j == i) {
                            modified_values_[k] = SCAL(1);  // Diagonal = 1
                        } else {
                            modified_values_[k] = SCAL(0);  // Off-diagonal = 0
                        }
                    } else if (!freedofs_->Test(j)) {
                        // Free row but column points to constrained DOF: zero it
                        modified_values_[k] = SCAL(0);
                    }
                }
            }

            sparsesolv::SparseMatrixView<SCAL> view(
                height_, width_,
                row_ptr_.data(),
                col_idx_.data(),
                modified_values_.data()
            );
            precond_->setup(view);
        } else {
            sparsesolv::SparseMatrixView<SCAL> view(
                height_, width_,
                row_ptr_.data(),
                col_idx_.data(),
                values.Data()
            );
            precond_->setup(view);
        }
    }

    /**
     * @brief Apply preconditioner: y = M^{-1} * x
     */
    void Mult(const BaseVector& x, BaseVector& y) const override {
        const SCAL* x_data = GetVectorData<SCAL>(x);
        SCAL* y_data = GetVectorData<SCAL>(y);

        if (freedofs_) {
            precond_->apply(x_data, y_data, height_);

            // Zero out constrained DOFs (OpenMP parallel)
            #pragma omp parallel for schedule(static)
            for (sparsesolv::index_t i = 0; i < height_; ++i) {
                if (!freedofs_->Test(i)) {
                    y_data[i] = SCAL(0);
                }
            }
        } else {
            precond_->apply(x_data, y_data, height_);
        }
    }

    void MultAdd(double s, const BaseVector& x, BaseVector& y) const override {
        const SCAL* x_data = GetVectorData<SCAL>(x);
        SCAL* y_data = GetVectorData<SCAL>(y);

        std::vector<SCAL> temp(height_);
        precond_->apply(x_data, temp.data(), height_);

        if (freedofs_) {
            #pragma omp parallel for schedule(static)
            for (sparsesolv::index_t i = 0; i < height_; ++i) {
                if (freedofs_->Test(i)) {
                    y_data[i] += s * temp[i];
                }
            }
        } else {
            #pragma omp parallel for schedule(static)
            for (sparsesolv::index_t i = 0; i < height_; ++i) {
                y_data[i] += s * temp[i];
            }
        }
    }

    int VHeight() const override { return static_cast<int>(height_); }
    int VWidth() const override { return static_cast<int>(width_); }

    AutoVector CreateRowVector() const override {
        return mat_->CreateRowVector();
    }

    AutoVector CreateColVector() const override {
        return mat_->CreateColVector();
    }

private:
    shared_ptr<SparseMatrix<SCAL>> mat_;
    shared_ptr<BitArray> freedofs_;
    sparsesolv::index_t height_;
    sparsesolv::index_t width_;
    shared_ptr<sparsesolv::SGSPreconditioner<SCAL>> precond_;

    std::vector<sparsesolv::index_t> row_ptr_;
    std::vector<sparsesolv::index_t> col_idx_;
    std::vector<SCAL> modified_values_;
};

// ============================================================================
// ILU Preconditioner with FreeDofs support
// ============================================================================

/**
 * @brief Incomplete LU Preconditioner for NGSolve
 *
 * Uses shifted ILU decomposition for improved stability.
 * Suitable for general (non-symmetric) matrices.
 * Supports Dirichlet boundary conditions through freedofs parameter.
 *
 * @tparam SCAL Scalar type (double or Complex)
 */
template<typename SCAL = double>
class SparseSolvILUPreconditioner : public BaseMatrix {
public:
    /**
     * @brief Construct ILU preconditioner from NGSolve sparse matrix
     * @param mat The sparse matrix to precondition
     * @param freedofs BitArray indicating free DOFs (nullptr = all DOFs are free)
     * @param shift Shift parameter for ILU decomposition (default: 1.05)
     */
    SparseSolvILUPreconditioner(shared_ptr<SparseMatrix<SCAL>> mat,
                                shared_ptr<BitArray> freedofs = nullptr,
                                double shift = 1.05)
        : mat_(mat)
        , freedofs_(freedofs ? make_shared<BitArray>(*freedofs) : nullptr)
        , shift_(shift)
        , height_(static_cast<sparsesolv::index_t>(mat->Height()))
        , width_(static_cast<sparsesolv::index_t>(mat->Width()))
        , precond_(std::make_shared<sparsesolv::ILUPreconditioner<SCAL>>(shift))
    {}

    /**
     * @brief Update the preconditioner (recompute factorization)
     */
    void Update() {
        const auto& firsti = mat_->GetFirstArray();
        const auto& colnr = mat_->GetColIndices();
        const auto& values = mat_->GetValues();

        const size_t nrows = firsti.Size();
        const size_t nnz = colnr.Size();

        // Convert row pointers (OpenMP parallel)
        row_ptr_.resize(nrows);
        #pragma omp parallel for schedule(static)
        for (ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(nrows); ++i) {
            row_ptr_[i] = static_cast<sparsesolv::index_t>(firsti[i]);
        }

        // Convert column indices (OpenMP parallel)
        col_idx_.resize(nnz);
        #pragma omp parallel for schedule(static)
        for (ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(nnz); ++i) {
            col_idx_[i] = static_cast<sparsesolv::index_t>(colnr[i]);
        }

        // Create modified values array if freedofs is specified
        // For constrained DOFs, set diagonal to 1 and off-diagonals to 0
        // Also zero out column entries pointing to constrained DOFs
        if (freedofs_) {
            modified_values_.resize(values.Size());

            // Copy values (OpenMP parallel)
            #pragma omp parallel for schedule(static)
            for (ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(values.Size()); ++i) {
                modified_values_[i] = values[i];
            }

            // Modify matrix entries for constrained DOFs (OpenMP parallel by row)
            #pragma omp parallel for schedule(dynamic, 64)
            for (sparsesolv::index_t i = 0; i < height_; ++i) {
                for (sparsesolv::index_t k = row_ptr_[i]; k < row_ptr_[i + 1]; ++k) {
                    sparsesolv::index_t j = col_idx_[k];
                    if (!freedofs_->Test(i)) {
                        // Constrained row: set to identity row
                        if (j == i) {
                            modified_values_[k] = SCAL(1);  // Diagonal = 1
                        } else {
                            modified_values_[k] = SCAL(0);  // Off-diagonal = 0
                        }
                    } else if (!freedofs_->Test(j)) {
                        // Free row but column points to constrained DOF: zero it
                        modified_values_[k] = SCAL(0);
                    }
                }
            }

            sparsesolv::SparseMatrixView<SCAL> view(
                height_, width_,
                row_ptr_.data(),
                col_idx_.data(),
                modified_values_.data()
            );
            precond_->setup(view);
        } else {
            sparsesolv::SparseMatrixView<SCAL> view(
                height_, width_,
                row_ptr_.data(),
                col_idx_.data(),
                values.Data()
            );
            precond_->setup(view);
        }
    }

    /**
     * @brief Apply preconditioner: y = M^{-1} * x
     */
    void Mult(const BaseVector& x, BaseVector& y) const override {
        const SCAL* x_data = GetVectorData<SCAL>(x);
        SCAL* y_data = GetVectorData<SCAL>(y);

        if (freedofs_) {
            // Apply preconditioner
            precond_->apply(x_data, y_data, height_);

            // Zero out constrained DOFs (OpenMP parallel)
            #pragma omp parallel for schedule(static)
            for (sparsesolv::index_t i = 0; i < height_; ++i) {
                if (!freedofs_->Test(i)) {
                    y_data[i] = SCAL(0);
                }
            }
        } else {
            precond_->apply(x_data, y_data, height_);
        }
    }

    void MultAdd(double s, const BaseVector& x, BaseVector& y) const override {
        const SCAL* x_data = GetVectorData<SCAL>(x);
        SCAL* y_data = GetVectorData<SCAL>(y);

        std::vector<SCAL> temp(height_);
        precond_->apply(x_data, temp.data(), height_);

        if (freedofs_) {
            #pragma omp parallel for schedule(static)
            for (sparsesolv::index_t i = 0; i < height_; ++i) {
                if (freedofs_->Test(i)) {
                    y_data[i] += s * temp[i];
                }
            }
        } else {
            #pragma omp parallel for schedule(static)
            for (sparsesolv::index_t i = 0; i < height_; ++i) {
                y_data[i] += s * temp[i];
            }
        }
    }

    int VHeight() const override { return static_cast<int>(height_); }
    int VWidth() const override { return static_cast<int>(width_); }

    AutoVector CreateRowVector() const override {
        return mat_->CreateRowVector();
    }

    AutoVector CreateColVector() const override {
        return mat_->CreateColVector();
    }

    double GetShift() const { return shift_; }

    void SetShift(double shift) {
        shift_ = shift;
        precond_->set_shift_parameter(shift);
    }

private:
    shared_ptr<SparseMatrix<SCAL>> mat_;
    shared_ptr<BitArray> freedofs_;
    double shift_;
    sparsesolv::index_t height_;
    sparsesolv::index_t width_;
    shared_ptr<sparsesolv::ILUPreconditioner<SCAL>> precond_;

    std::vector<sparsesolv::index_t> row_ptr_;
    std::vector<sparsesolv::index_t> col_idx_;
    std::vector<SCAL> modified_values_;
};

// Type aliases for convenience
using ICPreconditioner = SparseSolvICPreconditioner<double>;
using ILUPreconditioner = SparseSolvILUPreconditioner<double>;
using SGSPreconditioner = SparseSolvSGSPreconditioner<double>;

} // namespace ngla

#endif // NGSOLVE_SPARSESOLV_PRECOND_HPP
