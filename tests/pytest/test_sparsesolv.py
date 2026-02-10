"""
Tests for SparseSolv iterative solvers and preconditioners.

Tests SparseSolvSolver (ICCG, ICMRTR, SGSMRTR) and standalone
preconditioners (IC, ILU, SGS) for use with NGSolve's Krylov solvers.
"""

import pytest
from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from ngsolve import *
from ngsolve.la import ICPreconditioner, ILUPreconditioner, SGSPreconditioner
from ngsolve.la import SparseSolvSolver
from ngsolve.krylovspace import CGSolver


# ============================================================================
# Fixtures
# ============================================================================

@pytest.fixture
def poisson_2d():
    """2D Poisson problem with exact solution."""
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.1))
    fes = H1(mesh, order=2, dirichlet="bottom|right|top|left")
    u, v = fes.TnT()

    exact = sin(pi * x) * sin(pi * y)
    f_rhs = 2 * pi * pi * sin(pi * x) * sin(pi * y)

    a = BilinearForm(fes)
    a += grad(u) * grad(v) * dx
    a.Assemble()

    f = LinearForm(fes)
    f += f_rhs * v * dx
    f.Assemble()

    gfu_bc = GridFunction(fes)
    gfu_bc.Set(exact, BND)
    f.vec.data -= a.mat * gfu_bc.vec

    return mesh, fes, a, f, gfu_bc, exact


@pytest.fixture
def poisson_3d():
    """3D Poisson problem with exact solution."""
    mesh = Mesh(unit_cube.GenerateMesh(maxh=0.3))
    fes = H1(mesh, order=2, dirichlet="bottom|right|top|left|front|back")
    u, v = fes.TnT()

    exact = sin(pi * x) * sin(pi * y) * sin(pi * z)
    f_rhs = 3 * pi * pi * sin(pi * x) * sin(pi * y) * sin(pi * z)

    a = BilinearForm(fes)
    a += grad(u) * grad(v) * dx
    a.Assemble()

    f = LinearForm(fes)
    f += f_rhs * v * dx
    f.Assemble()

    gfu_bc = GridFunction(fes)
    gfu_bc.Set(exact, BND)
    f.vec.data -= a.mat * gfu_bc.vec

    return mesh, fes, a, f, gfu_bc, exact


# ============================================================================
# SparseSolvSolver tests
# ============================================================================

@pytest.mark.parametrize("method", ["ICCG", "ICMRTR", "SGSMRTR"])
def test_sparsesolv_solver_2d_poisson(poisson_2d, method):
    """All solver methods converge on 2D Poisson."""
    mesh, fes, a, f, gfu_bc, exact = poisson_2d

    solver = SparseSolvSolver(a.mat, method=method,
                               freedofs=fes.FreeDofs(),
                               tol=1e-10, maxiter=2000)

    gfu = GridFunction(fes)
    gfu.vec.data = solver * f.vec
    gfu.vec.data += gfu_bc.vec

    error = sqrt(Integrate((gfu - exact) ** 2, mesh))
    result = solver.last_result
    print(f"{method}: iterations={result.iterations}, L2 error={error:.2e}")
    assert result.converged
    assert result.iterations < 100
    assert error < 1e-3


def test_sparsesolv_solver_3d_poisson(poisson_3d):
    """ICCG converges on 3D Poisson."""
    mesh, fes, a, f, gfu_bc, exact = poisson_3d

    solver = SparseSolvSolver(a.mat, method="ICCG",
                               freedofs=fes.FreeDofs(),
                               tol=1e-8, maxiter=5000)

    gfu = GridFunction(fes)
    gfu.vec.data = solver * f.vec
    gfu.vec.data += gfu_bc.vec

    error = sqrt(Integrate((gfu - exact) ** 2, mesh))
    result = solver.last_result
    print(f"ICCG 3D: iterations={result.iterations}, L2 error={error:.2e}")
    assert result.converged
    assert error < 1e-2


def test_sparsesolv_solver_vs_direct(poisson_2d):
    """ICCG matches direct solver solution."""
    mesh, fes, a, f, gfu_bc, exact = poisson_2d

    gfu_direct = GridFunction(fes)
    gfu_direct.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec

    solver = SparseSolvSolver(a.mat, method="ICCG",
                               freedofs=fes.FreeDofs(),
                               tol=1e-12, maxiter=2000)
    gfu_iccg = GridFunction(fes)
    gfu_iccg.vec.data = solver * f.vec

    diff = gfu_iccg.vec.CreateVector()
    diff.data = gfu_iccg.vec - gfu_direct.vec
    rel_err = Norm(diff) / Norm(gfu_direct.vec)
    print(f"Relative error vs direct: {rel_err:.2e}")
    assert rel_err < 1e-6


# ============================================================================
# Preconditioner tests (with NGSolve CGSolver)
# ============================================================================

@pytest.mark.parametrize("PreClass,kwargs", [
    (ICPreconditioner, {"shift": 1.05}),
    (ILUPreconditioner, {"shift": 1.05}),
    (SGSPreconditioner, {}),
])
def test_preconditioners_with_ngsolve_cg(poisson_2d, PreClass, kwargs):
    """Preconditioners work with NGSolve's CGSolver."""
    mesh, fes, a, f, gfu_bc, exact = poisson_2d

    pre = PreClass(a.mat, freedofs=fes.FreeDofs(), **kwargs)
    pre.Update()

    inv = CGSolver(a.mat, pre, printrates=False, tol=1e-10, maxiter=2000)

    gfu = GridFunction(fes)
    gfu.vec.data = inv * f.vec
    gfu.vec.data += gfu_bc.vec

    error = sqrt(Integrate((gfu - exact) ** 2, mesh))
    print(f"{PreClass.__name__}: L2 error={error:.2e}")
    assert error < 1e-3


# ============================================================================
# Auto-shift and diagonal scaling (3D curl-curl)
# ============================================================================

def test_auto_shift_curl_curl():
    """Auto-shift + diagonal scaling for 3D curl-curl with nograds=True."""
    from netgen.occ import Box, Pnt

    box = Box(Pnt(0, 0, 0), Pnt(1, 1, 1))
    for face in box.faces:
        face.name = "outer"
    mesh = box.GenerateMesh(maxh=0.4)

    fes = HCurl(mesh, order=1, dirichlet="outer", nograds=True)
    u, v = fes.TnT()

    a = BilinearForm(fes)
    a += curl(u) * curl(v) * dx
    a.Assemble()

    f = LinearForm(fes)
    f += CF((0, 0, 1)) * v * dx
    f.Assemble()

    solver = SparseSolvSolver(a.mat, method="ICCG",
                               freedofs=fes.FreeDofs(),
                               tol=1e-8, maxiter=2000, shift=1.0)
    solver.auto_shift = True
    solver.diagonal_scaling = True

    gfu = GridFunction(fes)
    gfu.vec.data = solver * f.vec
    result = solver.last_result
    print(f"Auto-shift curl-curl: iterations={result.iterations}")
    assert result.converged
    assert result.iterations < 100


# ============================================================================
# Residual history and best-result tracking
# ============================================================================

def test_residual_history(poisson_2d):
    """Residual history is recorded when enabled."""
    _, fes, a, f, _, _ = poisson_2d

    solver = SparseSolvSolver(a.mat, method="ICCG",
                               freedofs=fes.FreeDofs(),
                               tol=1e-10, maxiter=2000,
                               save_residual_history=True)
    gfu = GridFunction(fes)
    result = solver.Solve(f.vec, gfu.vec)

    assert result.converged
    assert len(result.residual_history) > 0
    assert result.residual_history[-1] < result.residual_history[0]


def test_no_residual_history(poisson_2d):
    """Residual history is empty when disabled."""
    _, fes, a, f, _, _ = poisson_2d

    solver = SparseSolvSolver(a.mat, method="ICCG",
                               freedofs=fes.FreeDofs(),
                               save_residual_history=False)
    gfu = GridFunction(fes)
    result = solver.Solve(f.vec, gfu.vec)

    assert len(result.residual_history) == 0


# ============================================================================
# Property accessors
# ============================================================================

def test_properties(poisson_2d):
    """Property getters and setters work correctly."""
    _, fes, a, _, _, _ = poisson_2d

    solver = SparseSolvSolver(a.mat, method="ICCG",
                               freedofs=fes.FreeDofs())

    assert solver.method == "ICCG"
    assert abs(solver.tol - 1e-10) < 1e-15
    assert solver.maxiter == 1000
    assert abs(solver.shift - 1.05) < 1e-10
    assert solver.save_best_result is True
    assert solver.save_residual_history is False
    assert solver.auto_shift is False
    assert solver.diagonal_scaling is False

    solver.method = "SGSMRTR"
    assert solver.method == "SGSMRTR"

    solver.tol = 1e-8
    assert abs(solver.tol - 1e-8) < 1e-15

    solver.maxiter = 500
    assert solver.maxiter == 500

    solver.auto_shift = True
    assert solver.auto_shift is True

    solver.diagonal_scaling = True
    assert solver.diagonal_scaling is True


def test_invalid_method(poisson_2d):
    """Invalid method raises RuntimeError."""
    _, fes, a, f, _, _ = poisson_2d

    solver = SparseSolvSolver(a.mat, method="INVALID",
                               freedofs=fes.FreeDofs())
    gfu = GridFunction(fes)
    with pytest.raises(RuntimeError):
        solver.Solve(f.vec, gfu.vec)


# ============================================================================
# Operator interface
# ============================================================================

def test_operator_interface(poisson_2d):
    """solver * vec gives same result as solver.Solve()."""
    _, fes, a, f, _, _ = poisson_2d

    solver = SparseSolvSolver(a.mat, method="ICCG",
                               freedofs=fes.FreeDofs(),
                               tol=1e-10, maxiter=2000)

    gfu1 = GridFunction(fes)
    gfu1.vec.data = solver * f.vec

    gfu2 = GridFunction(fes)
    solver.Solve(f.vec, gfu2.vec)

    diff = gfu1.vec.CreateVector()
    diff.data = gfu1.vec - gfu2.vec
    assert Norm(diff) / max(Norm(gfu1.vec), 1e-30) < 1e-6


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
