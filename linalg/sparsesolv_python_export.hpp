/**
 * @file sparsesolv_python_export.hpp
 * @brief Pybind11 export functions for SparseSolv NGSolve integration
 *
 * Provides factory-based Python bindings following NGSolve's QMRSolver/CGSolver
 * pattern: typed classes with D/C suffix + factory functions with auto-dispatch
 * via mat->IsComplex().
 *
 * Usage in python_linalg.cpp:
 *   #include "sparsesolv_python_export.hpp"
 *
 *   // In ExportNgla():
 *   ExportSparseSolvBindings(m);
 *
 * Based on JP-MARs/SparseSolv (https://github.com/JP-MARs/SparseSolv)
 */

#ifndef NGSOLVE_SPARSESOLV_PYTHON_EXPORT_HPP
#define NGSOLVE_SPARSESOLV_PYTHON_EXPORT_HPP

#include "sparsesolv_precond.hpp"
#include <type_traits>

namespace ngla {

// ============================================================================
// Internal: SparseSolvResult (non-templated, called once)
// ============================================================================

inline void ExportSparseSolvResult_impl(py::module& m) {
  py::class_<SparseSolvResult>(m, "SparseSolvResult",
    R"raw_string(
Result of a SparseSolv iterative solve.

Attributes:

converged : bool
  Whether the solver converged within tolerance.

iterations : int
  Number of iterations performed.

final_residual : float
  Final relative residual. If save_best_result was enabled and the solver
  did not converge, this is the best residual achieved during iteration.

residual_history : list[float]
  Residual at each iteration (only if save_residual_history=True).

)raw_string")
    .def_readonly("converged", &SparseSolvResult::converged,
        "Whether the solver converged within tolerance")
    .def_readonly("iterations", &SparseSolvResult::iterations,
        "Number of iterations performed")
    .def_readonly("final_residual", &SparseSolvResult::final_residual,
        "Final relative residual (or best residual if save_best_result enabled)")
    .def_readonly("residual_history", &SparseSolvResult::residual_history,
        "Residual at each iteration (if save_residual_history enabled)")
    .def("__repr__", [](const SparseSolvResult& r) {
      return string("SparseSolvResult(converged=") +
             (r.converged ? "True" : "False") +
             ", iterations=" + std::to_string(r.iterations) +
             ", residual=" + std::to_string(r.final_residual) + ")";
    });
}

// ============================================================================
// Internal: Typed class registration (D/C suffix)
// ============================================================================

template<typename SCAL>
void ExportSparseSolvTyped(py::module& m, const std::string& suffix) {

  // IC Preconditioner (ICPreconditionerD / ICPreconditionerC)
  {
    std::string cls_name = "ICPreconditioner" + suffix;
    py::class_<SparseSolvICPreconditioner<SCAL>,
               shared_ptr<SparseSolvICPreconditioner<SCAL>>,
               BaseMatrix>(m, cls_name.c_str(),
               "Incomplete Cholesky preconditioner (typed). Use ICPreconditioner() factory instead.")
      .def(py::init([](shared_ptr<SparseMatrix<SCAL>> mat,
                       py::object freedofs, double shift) {
        shared_ptr<BitArray> sp_freedofs = nullptr;
        if (!freedofs.is_none())
          sp_freedofs = py::cast<shared_ptr<BitArray>>(freedofs);
        return make_shared<SparseSolvICPreconditioner<SCAL>>(mat, sp_freedofs, shift);
      }), py::arg("mat"), py::arg("freedofs") = py::none(), py::arg("shift") = 1.05)
      .def("Update", &SparseSolvICPreconditioner<SCAL>::Update,
          "Update preconditioner (recompute factorization after matrix change)")
      .def_property("shift",
          &SparseSolvICPreconditioner<SCAL>::GetShift,
          &SparseSolvICPreconditioner<SCAL>::SetShift,
          "Shift parameter for IC decomposition");
  }

  // SGS Preconditioner (SGSPreconditionerD / SGSPreconditionerC)
  {
    std::string cls_name = "SGSPreconditioner" + suffix;
    py::class_<SparseSolvSGSPreconditioner<SCAL>,
               shared_ptr<SparseSolvSGSPreconditioner<SCAL>>,
               BaseMatrix>(m, cls_name.c_str(),
               "Symmetric Gauss-Seidel preconditioner (typed). Use SGSPreconditioner() factory instead.")
      .def(py::init([](shared_ptr<SparseMatrix<SCAL>> mat,
                       py::object freedofs) {
        shared_ptr<BitArray> sp_freedofs = nullptr;
        if (!freedofs.is_none())
          sp_freedofs = py::cast<shared_ptr<BitArray>>(freedofs);
        return make_shared<SparseSolvSGSPreconditioner<SCAL>>(mat, sp_freedofs);
      }), py::arg("mat"), py::arg("freedofs") = py::none())
      .def("Update", &SparseSolvSGSPreconditioner<SCAL>::Update,
          "Update preconditioner (recompute after matrix change)");
  }

  // ILU Preconditioner (ILUPreconditionerD / ILUPreconditionerC)
  {
    std::string cls_name = "ILUPreconditioner" + suffix;
    py::class_<SparseSolvILUPreconditioner<SCAL>,
               shared_ptr<SparseSolvILUPreconditioner<SCAL>>,
               BaseMatrix>(m, cls_name.c_str(),
               "Incomplete LU preconditioner (typed). Use ILUPreconditioner() factory instead.")
      .def(py::init([](shared_ptr<SparseMatrix<SCAL>> mat,
                       py::object freedofs, double shift) {
        shared_ptr<BitArray> sp_freedofs = nullptr;
        if (!freedofs.is_none())
          sp_freedofs = py::cast<shared_ptr<BitArray>>(freedofs);
        return make_shared<SparseSolvILUPreconditioner<SCAL>>(mat, sp_freedofs, shift);
      }), py::arg("mat"), py::arg("freedofs") = py::none(), py::arg("shift") = 1.05)
      .def("Update", &SparseSolvILUPreconditioner<SCAL>::Update,
          "Update preconditioner (recompute factorization after matrix change)")
      .def_property("shift",
          &SparseSolvILUPreconditioner<SCAL>::GetShift,
          &SparseSolvILUPreconditioner<SCAL>::SetShift,
          "Shift parameter for ILU decomposition");
  }

  // SparseSolv Solver (SparseSolvSolverD / SparseSolvSolverC)
  {
    std::string cls_name = "SparseSolvSolver" + suffix;
    py::class_<SparseSolvSolver<SCAL>,
               shared_ptr<SparseSolvSolver<SCAL>>,
               BaseMatrix>(m, cls_name.c_str(),
               "SparseSolv iterative solver (typed). Use SparseSolvSolver() factory instead.")
      .def(py::init([](shared_ptr<SparseMatrix<SCAL>> mat,
                       const string& method, py::object freedofs,
                       double tol, int maxiter, double shift,
                       bool save_best_result, bool save_residual_history,
                       bool printrates) {
        shared_ptr<BitArray> sp_freedofs = nullptr;
        if (!freedofs.is_none())
          sp_freedofs = py::cast<shared_ptr<BitArray>>(freedofs);
        return make_shared<SparseSolvSolver<SCAL>>(
            mat, method, sp_freedofs, tol, maxiter, shift,
            save_best_result, save_residual_history, printrates);
      }),
          py::arg("mat"),
          py::arg("method") = "ICCG",
          py::arg("freedofs") = py::none(),
          py::arg("tol") = 1e-10,
          py::arg("maxiter") = 1000,
          py::arg("shift") = 1.05,
          py::arg("save_best_result") = true,
          py::arg("save_residual_history") = false,
          py::arg("printrates") = false)
      .def("Solve", [](SparseSolvSolver<SCAL>& self,
                       const BaseVector& rhs, BaseVector& sol) {
        return self.Solve(rhs, sol);
      }, py::arg("rhs"), py::arg("sol"),
         "Solve Ax = b. Returns SparseSolvResult with convergence info.")
      .def_property("method",
          &SparseSolvSolver<SCAL>::GetMethod,
          &SparseSolvSolver<SCAL>::SetMethod,
          "Solver method: ICCG, ICMRTR, SGSMRTR, CG, MRTR")
      .def_property("tol",
          &SparseSolvSolver<SCAL>::GetTolerance,
          &SparseSolvSolver<SCAL>::SetTolerance,
          "Relative convergence tolerance")
      .def_property("maxiter",
          &SparseSolvSolver<SCAL>::GetMaxIterations,
          &SparseSolvSolver<SCAL>::SetMaxIterations,
          "Maximum number of iterations")
      .def_property("shift",
          &SparseSolvSolver<SCAL>::GetShift,
          &SparseSolvSolver<SCAL>::SetShift,
          "Shift parameter for IC preconditioner")
      .def_property("save_best_result",
          &SparseSolvSolver<SCAL>::GetSaveBestResult,
          &SparseSolvSolver<SCAL>::SetSaveBestResult,
          "Track best solution during iteration")
      .def_property("save_residual_history",
          &SparseSolvSolver<SCAL>::GetSaveResidualHistory,
          &SparseSolvSolver<SCAL>::SetSaveResidualHistory,
          "Record residual at each iteration")
      .def_property("printrates",
          &SparseSolvSolver<SCAL>::GetPrintRates,
          &SparseSolvSolver<SCAL>::SetPrintRates,
          "Print convergence information after solve")
      .def_property("auto_shift",
          &SparseSolvSolver<SCAL>::GetAutoShift,
          &SparseSolvSolver<SCAL>::SetAutoShift,
          "Enable automatic shift adjustment for IC decomposition")
      .def_property("diagonal_scaling",
          &SparseSolvSolver<SCAL>::GetDiagonalScaling,
          &SparseSolvSolver<SCAL>::SetDiagonalScaling,
          "Enable diagonal scaling for IC preconditioner")
      .def_property("divergence_check",
          &SparseSolvSolver<SCAL>::GetDivergenceCheck,
          &SparseSolvSolver<SCAL>::SetDivergenceCheck,
          "Enable stagnation-based early termination")
      .def_property("divergence_threshold",
          &SparseSolvSolver<SCAL>::GetDivergenceThreshold,
          &SparseSolvSolver<SCAL>::SetDivergenceThreshold,
          "Multiplier for divergence detection (stop if residual > best * threshold)")
      .def_property("divergence_count",
          &SparseSolvSolver<SCAL>::GetDivergenceCount,
          &SparseSolvSolver<SCAL>::SetDivergenceCount,
          "Number of consecutive bad iterations before declaring divergence")
      .def_property_readonly("last_result",
          &SparseSolvSolver<SCAL>::GetLastResult,
          "Result from the last Solve() or Mult() call");
  }
}

// ============================================================================
// Internal: Helper to extract freedofs from py::object
// ============================================================================

inline shared_ptr<BitArray> ExtractFreeDofs(py::object freedofs) {
  if (freedofs.is_none()) return nullptr;
  return py::cast<shared_ptr<BitArray>>(freedofs);
}

// ============================================================================
// Internal: Factory functions with auto-dispatch via mat->IsComplex()
// ============================================================================

inline void ExportSparseSolvFactories(py::module& m) {

  // ---- ICPreconditioner factory ----
  m.def("ICPreconditioner", [](shared_ptr<BaseMatrix> mat,
                                py::object freedofs, double shift) {
    auto sp_freedofs = ExtractFreeDofs(freedofs);
    shared_ptr<BaseMatrix> result;
    if (mat->IsComplex()) {
      auto sp = dynamic_pointer_cast<SparseMatrix<Complex>>(mat);
      if (!sp) throw py::type_error("ICPreconditioner: expected SparseMatrix");
      auto p = make_shared<SparseSolvICPreconditioner<Complex>>(sp, sp_freedofs, shift);
      p->Update();
      result = p;
    } else {
      auto sp = dynamic_pointer_cast<SparseMatrix<double>>(mat);
      if (!sp) throw py::type_error("ICPreconditioner: expected SparseMatrix");
      auto p = make_shared<SparseSolvICPreconditioner<double>>(sp, sp_freedofs, shift);
      p->Update();
      result = p;
    }
    return result;
  },
  py::arg("mat"), py::arg("freedofs") = py::none(), py::arg("shift") = 1.05,
  R"raw_string(
Incomplete Cholesky (IC) Preconditioner.

Based on SparseSolv library by JP-MARs. Automatically detects real/complex
from the matrix type.

Example usage:

.. code-block:: python

    from ngsolve import *
    from ngsolve.krylovspace import CGSolver

    fes = H1(mesh, order=2, dirichlet="left|right|top|bottom")
    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += grad(u)*grad(v)*dx
    a.Assemble()

    pre = ICPreconditioner(a.mat, freedofs=fes.FreeDofs(), shift=1.1)
    inv = CGSolver(a.mat, pre, printrates=True, tol=1e-10)
    gfu.vec.data = inv * f.vec

Parameters:

mat : ngsolve.la.SparseMatrix
  The sparse matrix to precondition (must be SPD on free DOFs).
  Real or complex; type is auto-detected.

freedofs : ngsolve.BitArray, optional
  BitArray indicating free DOFs. Constrained DOFs are treated as identity.

shift : float
  Shift parameter for IC decomposition (default: 1.05).

)raw_string");

  // ---- SGSPreconditioner factory ----
  m.def("SGSPreconditioner", [](shared_ptr<BaseMatrix> mat,
                                  py::object freedofs) {
    auto sp_freedofs = ExtractFreeDofs(freedofs);
    shared_ptr<BaseMatrix> result;
    if (mat->IsComplex()) {
      auto sp = dynamic_pointer_cast<SparseMatrix<Complex>>(mat);
      if (!sp) throw py::type_error("SGSPreconditioner: expected SparseMatrix");
      auto p = make_shared<SparseSolvSGSPreconditioner<Complex>>(sp, sp_freedofs);
      p->Update();
      result = p;
    } else {
      auto sp = dynamic_pointer_cast<SparseMatrix<double>>(mat);
      if (!sp) throw py::type_error("SGSPreconditioner: expected SparseMatrix");
      auto p = make_shared<SparseSolvSGSPreconditioner<double>>(sp, sp_freedofs);
      p->Update();
      result = p;
    }
    return result;
  },
  py::arg("mat"), py::arg("freedofs") = py::none(),
  R"raw_string(
Symmetric Gauss-Seidel (SGS) Preconditioner.

Based on SparseSolv library by JP-MARs. Automatically detects real/complex
from the matrix type.

Example usage:

.. code-block:: python

    from ngsolve import *
    from ngsolve.krylovspace import CGSolver

    fes = H1(mesh, order=2, dirichlet="left|right|top|bottom")
    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += grad(u)*grad(v)*dx
    a.Assemble()

    pre = SGSPreconditioner(a.mat, freedofs=fes.FreeDofs())
    inv = CGSolver(a.mat, pre, printrates=True)
    gfu.vec.data = inv * f.vec

Parameters:

mat : ngsolve.la.SparseMatrix
  The sparse matrix to precondition (must be SPD on free DOFs).
  Real or complex; type is auto-detected.

freedofs : ngsolve.BitArray, optional
  BitArray indicating free DOFs. Constrained DOFs are treated as identity.

)raw_string");

  // ---- ILUPreconditioner factory ----
  m.def("ILUPreconditioner", [](shared_ptr<BaseMatrix> mat,
                                  py::object freedofs, double shift) {
    auto sp_freedofs = ExtractFreeDofs(freedofs);
    shared_ptr<BaseMatrix> result;
    if (mat->IsComplex()) {
      auto sp = dynamic_pointer_cast<SparseMatrix<Complex>>(mat);
      if (!sp) throw py::type_error("ILUPreconditioner: expected SparseMatrix");
      auto p = make_shared<SparseSolvILUPreconditioner<Complex>>(sp, sp_freedofs, shift);
      p->Update();
      result = p;
    } else {
      auto sp = dynamic_pointer_cast<SparseMatrix<double>>(mat);
      if (!sp) throw py::type_error("ILUPreconditioner: expected SparseMatrix");
      auto p = make_shared<SparseSolvILUPreconditioner<double>>(sp, sp_freedofs, shift);
      p->Update();
      result = p;
    }
    return result;
  },
  py::arg("mat"), py::arg("freedofs") = py::none(), py::arg("shift") = 1.05,
  R"raw_string(
Incomplete LU (ILU) Preconditioner.

Based on SparseSolv library by JP-MARs. Automatically detects real/complex
from the matrix type.

Suitable for general (non-symmetric) matrices. For symmetric positive definite
matrices, ICPreconditioner is more efficient.

Example usage:

.. code-block:: python

    from ngsolve import *
    from ngsolve.krylovspace import GMRESSolver

    fes = H1(mesh, order=2, dirichlet="left|right|top|bottom")
    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += grad(u)*grad(v)*dx + u*v*dx
    a.Assemble()

    pre = ILUPreconditioner(a.mat, freedofs=fes.FreeDofs(), shift=1.1)
    inv = GMRESSolver(a.mat, pre, printrates=True, tol=1e-10)
    gfu.vec.data = inv * f.vec

Parameters:

mat : ngsolve.la.SparseMatrix
  The sparse matrix to precondition.
  Real or complex; type is auto-detected.

freedofs : ngsolve.BitArray, optional
  BitArray indicating free DOFs. Constrained DOFs are treated as identity.

shift : float
  Shift parameter for ILU decomposition (default: 1.05).

)raw_string");

  // ---- SparseSolvSolver factory ----
  m.def("SparseSolvSolver", [](shared_ptr<BaseMatrix> mat,
                                 const string& method, py::object freedofs,
                                 double tol, int maxiter, double shift,
                                 bool save_best_result, bool save_residual_history,
                                 bool printrates) {
    auto sp_freedofs = ExtractFreeDofs(freedofs);
    shared_ptr<BaseMatrix> result;
    if (mat->IsComplex()) {
      auto sp = dynamic_pointer_cast<SparseMatrix<Complex>>(mat);
      if (!sp) throw py::type_error("SparseSolvSolver: expected SparseMatrix");
      result = make_shared<SparseSolvSolver<Complex>>(
          sp, method, sp_freedofs, tol, maxiter, shift,
          save_best_result, save_residual_history, printrates);
    } else {
      auto sp = dynamic_pointer_cast<SparseMatrix<double>>(mat);
      if (!sp) throw py::type_error("SparseSolvSolver: expected SparseMatrix");
      result = make_shared<SparseSolvSolver<double>>(
          sp, method, sp_freedofs, tol, maxiter, shift,
          save_best_result, save_residual_history, printrates);
    }
    return result;
  },
  py::arg("mat"),
  py::arg("method") = "ICCG",
  py::arg("freedofs") = py::none(),
  py::arg("tol") = 1e-10,
  py::arg("maxiter") = 1000,
  py::arg("shift") = 1.05,
  py::arg("save_best_result") = true,
  py::arg("save_residual_history") = false,
  py::arg("printrates") = false,
  R"raw_string(
Iterative solver using the SparseSolv library by JP-MARs.

Automatically detects real/complex from the matrix type.

Supports multiple solver methods for symmetric positive definite systems:
- ICCG: Conjugate Gradient + Incomplete Cholesky preconditioner
- ICMRTR: MRTR (Modified Residual Tri-diagonal Reduction) + IC preconditioner
- SGSMRTR: MRTR with built-in Symmetric Gauss-Seidel (split formula)
- CG: Conjugate Gradient without preconditioner
- MRTR: MRTR without preconditioner

Key features:
- save_best_result (default: True): tracks best solution during iteration.
  If the solver doesn't converge, the best solution found is returned.
- save_residual_history: records residual at every iteration for analysis.
- divergence_check: early termination when residual stagnates relative to best.
- FreeDofs support for Dirichlet boundary conditions.

Can be used as an inverse operator (BaseMatrix) or with Solve() method.

Example usage as inverse operator:

.. code-block:: python

    from ngsolve import *

    fes = H1(mesh, order=2, dirichlet="left|right|top|bottom")
    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += grad(u)*grad(v)*dx
    a.Assemble()
    f = LinearForm(fes)
    f += 1*v*dx
    f.Assemble()

    gfu = GridFunction(fes)
    solver = SparseSolvSolver(a.mat, method="ICCG",
                              freedofs=fes.FreeDofs(), tol=1e-10)
    gfu.vec.data = solver * f.vec

Example usage with Solve() for detailed results:

.. code-block:: python

    solver = SparseSolvSolver(a.mat, method="ICCG",
                              freedofs=fes.FreeDofs(), tol=1e-10,
                              save_residual_history=True)
    result = solver.Solve(f.vec, gfu.vec)
    print(f"Converged: {result.converged}")
    print(f"Iterations: {result.iterations}")

Parameters:

mat : ngsolve.la.SparseMatrix
  The sparse system matrix (must be SPD for CG/MRTR methods).
  Real or complex; type is auto-detected.

method : str
  Solver method. One of: "ICCG", "ICMRTR", "SGSMRTR", "CG", "MRTR".

freedofs : ngsolve.BitArray, optional
  BitArray indicating free DOFs. Constrained DOFs are treated as identity.

tol : float
  Relative convergence tolerance (default: 1e-10).

maxiter : int
  Maximum number of iterations (default: 1000).

shift : float
  Shift parameter for IC preconditioner (default: 1.05).

save_best_result : bool
  Track best solution during iteration (default: True).

save_residual_history : bool
  Record residual at each iteration (default: False).

printrates : bool
  Print convergence information after solve (default: False).

Properties (set after construction):

divergence_check : bool
  Enable stagnation-based early termination (default: False).
  When enabled, the solver stops if the residual remains worse than
  best_residual * divergence_threshold for divergence_count consecutive
  iterations.

divergence_threshold : float
  Multiplier for divergence detection (default: 1000.0).
  The solver counts an iteration as "bad" if
  residual >= best_residual * divergence_threshold.

divergence_count : int
  Number of consecutive bad iterations before declaring divergence
  (default: 100).

auto_shift : bool
  Enable automatic shift adjustment for IC decomposition (default: False).

diagonal_scaling : bool
  Enable diagonal scaling for IC preconditioner (default: False).

)raw_string");
}

// ============================================================================
// Public API: Single entry point for NGSolve integration
// ============================================================================

/**
 * @brief Register all SparseSolv Python bindings
 *
 * Call once from ExportNgla() in python_linalg.cpp.
 * Registers typed classes (ICPreconditionerD/C, etc.) and
 * factory functions (ICPreconditioner, etc.) with auto-dispatch.
 */
inline void ExportSparseSolvBindings(py::module& m) {
  ExportSparseSolvResult_impl(m);
  ExportSparseSolvTyped<double>(m, "D");
  ExportSparseSolvTyped<Complex>(m, "C");
  ExportSparseSolvFactories(m);
}

} // namespace ngla

#endif // NGSOLVE_SPARSESOLV_PYTHON_EXPORT_HPP
