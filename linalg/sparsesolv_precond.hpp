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
// Preconditioner Base Class
// ============================================================================

/**
 * @brief Base class for SparseSolv preconditioners in NGSolve
 *
 * Provides common functionality for all SparseSolv preconditioner wrappers:
 * - NGSolve SparseMatrix to SparseSolv SparseMatrixView conversion
 * - FreeDofs handling (constrained DOFs become identity rows/columns)
 * - Mult/MultAdd with FreeDofs zero-out
 * - VHeight/VWidth/CreateRowVector/CreateColVector
 *
 * Derived classes only need to implement:
 * - Constructor (create the specific sparsesolv preconditioner)
 * - Update() (call prepare_matrix_view() + precond setup)
 * - apply_precond() (delegate to the specific preconditioner)
 *
 * @tparam SCAL Scalar type (double or Complex)
 */
template<typename SCAL>
class SparseSolvPrecondBase : public BaseMatrix {
protected:
    shared_ptr<SparseMatrix<SCAL>> mat_;
    shared_ptr<BitArray> freedofs_;
    sparsesolv::index_t height_;
    sparsesolv::index_t width_;

    // CSR conversion buffers
    std::vector<sparsesolv::index_t> row_ptr_;
    std::vector<sparsesolv::index_t> col_idx_;
    std::vector<SCAL> modified_values_;

    SparseSolvPrecondBase(shared_ptr<SparseMatrix<SCAL>> mat,
                          shared_ptr<BitArray> freedofs)
        : mat_(mat)
        , freedofs_(freedofs ? make_shared<BitArray>(*freedofs) : nullptr)
        , height_(static_cast<sparsesolv::index_t>(mat->Height()))
        , width_(static_cast<sparsesolv::index_t>(mat->Width()))
    {}

    /**
     * @brief Convert NGSolve SparseMatrix to SparseSolv SparseMatrixView
     *
     * Handles FreeDofs: constrained DOFs get identity rows/columns
     * (diagonal = 1, off-diagonals = 0).
     */
    sparsesolv::SparseMatrixView<SCAL> prepare_matrix_view() {
        const auto& firsti = mat_->GetFirstArray();
        const auto& colnr = mat_->GetColIndices();
        const auto& values = mat_->GetValues();

        const size_t nrows = firsti.Size();
        const size_t nnz = colnr.Size();

        // Convert row pointers
        row_ptr_.resize(nrows);
        sparsesolv::parallel_for(static_cast<sparsesolv::index_t>(nrows), [&](sparsesolv::index_t i) {
            row_ptr_[i] = static_cast<sparsesolv::index_t>(firsti[i]);
        });

        // Convert column indices
        col_idx_.resize(nnz);
        sparsesolv::parallel_for(static_cast<sparsesolv::index_t>(nnz), [&](sparsesolv::index_t i) {
            col_idx_[i] = static_cast<sparsesolv::index_t>(colnr[i]);
        });

        if (freedofs_) {
            modified_values_.resize(values.Size());

            sparsesolv::parallel_for(static_cast<sparsesolv::index_t>(values.Size()), [&](sparsesolv::index_t i) {
                modified_values_[i] = values[i];
            });

            sparsesolv::parallel_for(height_, [&](sparsesolv::index_t i) {
                for (sparsesolv::index_t k = row_ptr_[i]; k < row_ptr_[i + 1]; ++k) {
                    sparsesolv::index_t j = col_idx_[k];
                    if (!freedofs_->Test(i)) {
                        modified_values_[k] = (j == i) ? SCAL(1) : SCAL(0);
                    } else if (!freedofs_->Test(j)) {
                        modified_values_[k] = SCAL(0);
                    }
                }
            });

            return sparsesolv::SparseMatrixView<SCAL>(
                height_, width_,
                row_ptr_.data(), col_idx_.data(), modified_values_.data()
            );
        } else {
            return sparsesolv::SparseMatrixView<SCAL>(
                height_, width_,
                row_ptr_.data(), col_idx_.data(), values.Data()
            );
        }
    }

    /**
     * @brief Apply the underlying SparseSolv preconditioner
     *
     * Derived classes implement this to delegate to their specific preconditioner.
     */
    virtual void apply_precond(const SCAL* x, SCAL* y) const = 0;

public:
    void Mult(const BaseVector& x, BaseVector& y) const override {
        const SCAL* x_data = GetVectorData<SCAL>(x);
        SCAL* y_data = GetVectorData<SCAL>(y);

        apply_precond(x_data, y_data);

        if (freedofs_) {
            sparsesolv::parallel_for(height_, [&](sparsesolv::index_t i) {
                if (!freedofs_->Test(i)) {
                    y_data[i] = SCAL(0);
                }
            });
        }
    }

    void MultAdd(double s, const BaseVector& x, BaseVector& y) const override {
        const SCAL* x_data = GetVectorData<SCAL>(x);
        SCAL* y_data = GetVectorData<SCAL>(y);

        std::vector<SCAL> temp(height_);
        apply_precond(x_data, temp.data());

        if (freedofs_) {
            sparsesolv::parallel_for(height_, [&](sparsesolv::index_t i) {
                if (freedofs_->Test(i)) {
                    y_data[i] += s * temp[i];
                }
            });
        } else {
            sparsesolv::parallel_for(height_, [&](sparsesolv::index_t i) {
                y_data[i] += s * temp[i];
            });
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
};

// ============================================================================
// IC Preconditioner
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
class SparseSolvICPreconditioner : public SparseSolvPrecondBase<SCAL> {
public:
    SparseSolvICPreconditioner(shared_ptr<SparseMatrix<SCAL>> mat,
                                shared_ptr<BitArray> freedofs = nullptr,
                                double shift = 1.05)
        : SparseSolvPrecondBase<SCAL>(mat, freedofs)
        , shift_(shift)
        , precond_(std::make_shared<sparsesolv::ICPreconditioner<SCAL>>(shift))
    {}

    void Update() {
        auto view = this->prepare_matrix_view();
        precond_->setup(view);
    }

    double GetShift() const { return shift_; }

    void SetShift(double shift) {
        shift_ = shift;
        precond_->set_shift_parameter(shift);
    }

protected:
    void apply_precond(const SCAL* x, SCAL* y) const override {
        precond_->apply(x, y, this->height_);
    }

private:
    double shift_;
    shared_ptr<sparsesolv::ICPreconditioner<SCAL>> precond_;
};

// ============================================================================
// SGS Preconditioner
// ============================================================================

/**
 * @brief Symmetric Gauss-Seidel Preconditioner for NGSolve
 *
 * Supports Dirichlet boundary conditions through freedofs parameter.
 *
 * @tparam SCAL Scalar type (double or Complex)
 */
template<typename SCAL = double>
class SparseSolvSGSPreconditioner : public SparseSolvPrecondBase<SCAL> {
public:
    SparseSolvSGSPreconditioner(shared_ptr<SparseMatrix<SCAL>> mat,
                                 shared_ptr<BitArray> freedofs = nullptr)
        : SparseSolvPrecondBase<SCAL>(mat, freedofs)
        , precond_(std::make_shared<sparsesolv::SGSPreconditioner<SCAL>>())
    {}

    void Update() {
        auto view = this->prepare_matrix_view();
        precond_->setup(view);
    }

protected:
    void apply_precond(const SCAL* x, SCAL* y) const override {
        precond_->apply(x, y, this->height_);
    }

private:
    shared_ptr<sparsesolv::SGSPreconditioner<SCAL>> precond_;
};

// ============================================================================
// ILU Preconditioner
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
class SparseSolvILUPreconditioner : public SparseSolvPrecondBase<SCAL> {
public:
    SparseSolvILUPreconditioner(shared_ptr<SparseMatrix<SCAL>> mat,
                                shared_ptr<BitArray> freedofs = nullptr,
                                double shift = 1.05)
        : SparseSolvPrecondBase<SCAL>(mat, freedofs)
        , shift_(shift)
        , precond_(std::make_shared<sparsesolv::ILUPreconditioner<SCAL>>(shift))
    {}

    void Update() {
        auto view = this->prepare_matrix_view();
        precond_->setup(view);
    }

    double GetShift() const { return shift_; }

    void SetShift(double shift) {
        shift_ = shift;
        precond_->set_shift_parameter(shift);
    }

protected:
    void apply_precond(const SCAL* x, SCAL* y) const override {
        precond_->apply(x, y, this->height_);
    }

private:
    double shift_;
    shared_ptr<sparsesolv::ILUPreconditioner<SCAL>> precond_;
};

// Type aliases for convenience
using ICPreconditioner = SparseSolvICPreconditioner<double>;
using ILUPreconditioner = SparseSolvILUPreconditioner<double>;
using SGSPreconditioner = SparseSolvSGSPreconditioner<double>;

// ============================================================================
// Solver Result (returned to Python)
// ============================================================================

/**
 * @brief Result of a SparseSolv iterative solve
 */
struct SparseSolvResult {
    bool converged = false;
    int iterations = 0;
    double final_residual = 0.0;
    std::vector<double> residual_history;
};

// ============================================================================
// SparseSolv Iterative Solver for NGSolve
// ============================================================================

/**
 * @brief Unified iterative solver using SparseSolv library
 *
 * Supports multiple solver methods:
 * - "ICCG": Conjugate Gradient with Incomplete Cholesky preconditioner
 * - "ICMRTR": MRTR with Incomplete Cholesky preconditioner
 * - "SGSMRTR": MRTR with built-in Symmetric Gauss-Seidel (split formula)
 * - "CG": Conjugate Gradient without preconditioner
 * - "MRTR": MRTR without preconditioner
 *
 * Key features:
 * - save_best_result: tracks best solution during iteration (default: true)
 *   If the solver doesn't converge, the best solution found is returned.
 * - save_residual_history: records residual at each iteration
 * - FreeDofs support for Dirichlet boundary conditions
 *
 * Can be used as an inverse operator (BaseMatrix) or with Solve() for
 * detailed result access.
 *
 * @code
 * from ngsolve import *
 *
 * # As inverse operator (like NGSolve's CGSolver)
 * solver = SparseSolvSolver(a.mat, method="ICCG",
 *                           freedofs=fes.FreeDofs(), tol=1e-10)
 * gfu.vec.data = solver * f.vec
 *
 * # Or with detailed results
 * result = solver.Solve(f.vec, gfu.vec)
 * print(f"Converged: {result.converged}, iterations: {result.iterations}")
 * @endcode
 *
 * @tparam SCAL Scalar type (double or Complex)
 */
template<typename SCAL = double>
class SparseSolvSolver : public BaseMatrix {
public:
    /**
     * @brief Construct solver
     * @param mat The sparse matrix (system matrix A)
     * @param method Solver method: "ICCG", "ICMRTR", "SGSMRTR", "CG", "MRTR"
     * @param freedofs BitArray indicating free DOFs (nullptr = all free)
     * @param tol Relative tolerance for convergence (default: 1e-10)
     * @param maxiter Maximum number of iterations (default: 1000)
     * @param shift Shift parameter for IC/ILU preconditioner (default: 1.05)
     * @param save_best_result Save best solution during iteration (default: true)
     * @param save_residual_history Save residual at each iteration (default: false)
     * @param printrates Print convergence information (default: false)
     */
    SparseSolvSolver(shared_ptr<SparseMatrix<SCAL>> mat,
                      const string& method = "ICCG",
                      shared_ptr<BitArray> freedofs = nullptr,
                      double tol = 1e-10,
                      int maxiter = 1000,
                      double shift = 1.05,
                      bool save_best_result = true,
                      bool save_residual_history = false,
                      bool printrates = false)
        : mat_(mat)
        , method_(method)
        , freedofs_(freedofs ? make_shared<BitArray>(*freedofs) : nullptr)
        , height_(static_cast<sparsesolv::index_t>(mat->Height()))
        , width_(static_cast<sparsesolv::index_t>(mat->Width()))
        , printrates_(printrates)
    {
        config_.tolerance = tol;
        config_.max_iterations = maxiter;
        config_.shift_parameter = shift;
        config_.save_best_result = save_best_result;
        config_.save_residual_history = save_residual_history;
    }

    /**
     * @brief Solve Ax = b, x initialized to zero (BaseMatrix interface)
     *
     * Use as: gfu.vec.data = solver * f.vec
     */
    void Mult(const BaseVector& x, BaseVector& y) const override {
        y = 0.0;
        Solve(x, y);
    }

    void MultAdd(double s, const BaseVector& x, BaseVector& y) const override {
        auto temp = y.CreateVector();
        temp = 0.0;
        Solve(x, *temp);
        y += s * *temp;
    }

    /**
     * @brief Solve Ax = b with initial guess, returning detailed results
     *
     * @param rhs Right-hand side vector b
     * @param sol Solution vector x (initial guess on input, solution on output)
     * @return SparseSolvResult with convergence info, iterations, residual history
     */
    SparseSolvResult Solve(const BaseVector& rhs, BaseVector& sol) const {
        // Prepare matrix view (with FreeDofs handling)
        auto view = prepare_matrix();

        // Get vector data
        const SCAL* b = GetVectorData<SCAL>(rhs);
        SCAL* x = GetVectorData<SCAL>(sol);

        // For FreeDofs: work on copies with constrained DOFs zeroed
        std::vector<SCAL> b_mod;
        std::vector<SCAL> x_mod;
        const SCAL* b_ptr = b;
        SCAL* x_ptr = x;

        if (freedofs_) {
            b_mod.resize(height_);
            x_mod.resize(height_);
            sparsesolv::parallel_for(height_, [&](sparsesolv::index_t i) {
                if (freedofs_->Test(i)) {
                    b_mod[i] = b[i];
                    x_mod[i] = x[i];
                } else {
                    b_mod[i] = SCAL(0);
                    x_mod[i] = SCAL(0);
                }
            });
            b_ptr = b_mod.data();
            x_ptr = x_mod.data();
        }

        // Configure solver
        sparsesolv::SolverConfig config = config_;
        if (printrates_) {
            config.save_residual_history = true;
        }

        // Dispatch to appropriate solver
        sparsesolv::SolverResult result;

        if (method_ == "ICCG" || method_ == "iccg") {
            result = sparsesolv::solve_iccg(view, b_ptr, x_ptr, height_, config);
        } else if (method_ == "ICMRTR" || method_ == "icmrtr") {
            result = sparsesolv::solve_icmrtr(view, b_ptr, x_ptr, height_, config);
        } else if (method_ == "SGSMRTR" || method_ == "sgsmrtr") {
            result = sparsesolv::solve_sgsmrtr(view, b_ptr, x_ptr, height_, config);
        } else if (method_ == "CG" || method_ == "cg") {
            sparsesolv::CGSolver<SCAL> solver;
            solver.set_config(config);
            result = solver.solve(view, b_ptr, x_ptr, height_, nullptr);
        } else if (method_ == "MRTR" || method_ == "mrtr") {
            sparsesolv::MRTRSolver<SCAL> solver;
            solver.set_config(config);
            result = solver.solve(view, b_ptr, x_ptr, height_, nullptr);
        } else {
            throw std::runtime_error(
                "SparseSolvSolver: Unknown method '" + method_ +
                "'. Available: ICCG, ICMRTR, SGSMRTR, CG, MRTR");
        }

        // Copy solution back (only free DOFs)
        if (freedofs_) {
            sparsesolv::parallel_for(height_, [&](sparsesolv::index_t i) {
                if (freedofs_->Test(i)) {
                    x[i] = x_ptr[i];
                }
            });
        }

        // Print convergence info if requested
        if (printrates_ && !result.residual_history.empty()) {
            for (size_t i = 0; i < result.residual_history.size(); ++i) {
                std::cout << method_ << " iteration " << i
                          << ", residual = " << result.residual_history[i]
                          << std::endl;
            }
            if (result.converged) {
                std::cout << method_ << " converged in "
                          << result.iterations << " iterations, residual = "
                          << result.final_residual << std::endl;
            } else {
                std::cout << method_ << " NOT converged after "
                          << result.iterations << " iterations, residual = "
                          << result.final_residual << std::endl;
            }
        }

        // Build and store result
        last_result_.converged = result.converged;
        last_result_.iterations = result.iterations;
        last_result_.final_residual = result.final_residual;
        last_result_.residual_history = std::move(result.residual_history);

        return last_result_;
    }

    int VHeight() const override { return static_cast<int>(height_); }
    int VWidth() const override { return static_cast<int>(width_); }
    AutoVector CreateRowVector() const override { return mat_->CreateRowVector(); }
    AutoVector CreateColVector() const override { return mat_->CreateColVector(); }

    // Property accessors
    const string& GetMethod() const { return method_; }
    void SetMethod(const string& method) { method_ = method; }
    double GetTolerance() const { return config_.tolerance; }
    void SetTolerance(double tol) { config_.tolerance = tol; }
    int GetMaxIterations() const { return config_.max_iterations; }
    void SetMaxIterations(int maxiter) { config_.max_iterations = maxiter; }
    double GetShift() const { return config_.shift_parameter; }
    void SetShift(double shift) { config_.shift_parameter = shift; }
    bool GetSaveBestResult() const { return config_.save_best_result; }
    void SetSaveBestResult(bool save) { config_.save_best_result = save; }
    bool GetSaveResidualHistory() const { return config_.save_residual_history; }
    void SetSaveResidualHistory(bool save) { config_.save_residual_history = save; }
    bool GetPrintRates() const { return printrates_; }
    void SetPrintRates(bool print) { printrates_ = print; }
    bool GetAutoShift() const { return config_.auto_shift; }
    void SetAutoShift(bool enable) { config_.auto_shift = enable; }
    bool GetDiagonalScaling() const { return config_.diagonal_scaling; }
    void SetDiagonalScaling(bool enable) { config_.diagonal_scaling = enable; }

    // Divergence detection
    bool GetDivergenceCheck() const {
        return config_.divergence_check == sparsesolv::DivergenceCheck::StagnationCount;
    }
    void SetDivergenceCheck(bool enable) {
        config_.divergence_check = enable
            ? sparsesolv::DivergenceCheck::StagnationCount
            : sparsesolv::DivergenceCheck::None;
    }
    double GetDivergenceThreshold() const { return config_.divergence_threshold; }
    void SetDivergenceThreshold(double threshold) { config_.divergence_threshold = threshold; }
    int GetDivergenceCount() const { return config_.divergence_count; }
    void SetDivergenceCount(int count) { config_.divergence_count = count; }

    const SparseSolvResult& GetLastResult() const { return last_result_; }

private:
    /**
     * @brief Prepare SparseMatrixView with FreeDofs handling
     *
     * Converts NGSolve's SparseMatrix to SparseSolv's SparseMatrixView.
     * For constrained DOFs: diagonal = 1, off-diagonals = 0 (identity row/col).
     */
    sparsesolv::SparseMatrixView<SCAL> prepare_matrix() const {
        const auto& firsti = mat_->GetFirstArray();
        const auto& colnr = mat_->GetColIndices();
        const auto& values = mat_->GetValues();

        const size_t nrows = firsti.Size();
        const size_t nnz = colnr.Size();

        // Convert row pointers
        row_ptr_.resize(nrows);
        sparsesolv::parallel_for(static_cast<sparsesolv::index_t>(nrows), [&](sparsesolv::index_t i) {
            row_ptr_[i] = static_cast<sparsesolv::index_t>(firsti[i]);
        });

        // Convert column indices
        col_idx_.resize(nnz);
        sparsesolv::parallel_for(static_cast<sparsesolv::index_t>(nnz), [&](sparsesolv::index_t i) {
            col_idx_[i] = static_cast<sparsesolv::index_t>(colnr[i]);
        });

        if (freedofs_) {
            modified_values_.resize(values.Size());

            sparsesolv::parallel_for(static_cast<sparsesolv::index_t>(values.Size()), [&](sparsesolv::index_t i) {
                modified_values_[i] = values[i];
            });

            sparsesolv::parallel_for(height_, [&](sparsesolv::index_t i) {
                for (sparsesolv::index_t k = row_ptr_[i]; k < row_ptr_[i + 1]; ++k) {
                    sparsesolv::index_t j = col_idx_[k];
                    if (!freedofs_->Test(i)) {
                        modified_values_[k] = (j == i) ? SCAL(1) : SCAL(0);
                    } else if (!freedofs_->Test(j)) {
                        modified_values_[k] = SCAL(0);
                    }
                }
            });

            return sparsesolv::SparseMatrixView<SCAL>(
                height_, width_,
                row_ptr_.data(), col_idx_.data(), modified_values_.data()
            );
        } else {
            return sparsesolv::SparseMatrixView<SCAL>(
                height_, width_,
                row_ptr_.data(), col_idx_.data(), values.Data()
            );
        }
    }

    shared_ptr<SparseMatrix<SCAL>> mat_;
    string method_;
    shared_ptr<BitArray> freedofs_;
    sparsesolv::SolverConfig config_;
    sparsesolv::index_t height_;
    sparsesolv::index_t width_;
    bool printrates_;

    // Mutable for use in const Mult/Solve methods
    mutable std::vector<sparsesolv::index_t> row_ptr_;
    mutable std::vector<sparsesolv::index_t> col_idx_;
    mutable std::vector<SCAL> modified_values_;
    mutable SparseSolvResult last_result_;
};

} // namespace ngla

#endif // NGSOLVE_SPARSESOLV_PRECOND_HPP
