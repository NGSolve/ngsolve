import pytest
from ngsolve import *
from ngsolve.comp import IntegrationRuleSpace
from netgen.geom2d import unit_square
import numpy as np


def mk_2d_mesh():
    return Mesh(unit_square.GenerateMesh())

@pytest.fixture
def mesh2d():
    return mk_2d_mesh()


def _test_interpolation(space, target_expr, check):
    gf = GridFunction(space)
    gf.Interpolate(target_expr)
    assert check(gf)
    gf.Set(target_expr)
    assert check(gf)
    gf.Interpolate(gf)
    assert check(gf)


def test_IRSpace2d(mesh2d):
    mesh = mesh2d
    space = IntegrationRuleSpace(mesh, order=2)

    def _expr(x, y):
        return 6 * x**2 + 4 * y**2

    eval_points = GridFunction(space**2)
    eval_points.Interpolate((x, y))
    q_points = eval_points.vec.FV().NumPy()
    q_points = q_points.reshape(2, (int(q_points.shape[0] / 2)))

    def _check(gf):
        return np.allclose(gf.vec.FV().NumPy(), _expr(q_points[0], q_points[1]), atol=1e-10, rtol=1e-10)

    _test_interpolation(space, _expr(x, y), _check)


def test_MatrixValuedL2(mesh2d):
    mesh = mesh2d
    space = MatrixValued(L2(mesh, order=2), dim=2)

    def _expr(x, y):
        return ((6 * x**2 + 4 * y**2, y**1 - 1.5 * x**2),
                (2 * x**2 + 5 * y**2, y**2 - 1.5 * x**1))

    def _check(gf):
        return np.allclose(gf(mesh(0.5, 0.5)), np.array(_expr(0.5, 0.5)).flatten(), atol=1e-10, rtol=1e-10)

    _test_interpolation(space, _expr(x, y), _check)


def test_SymMatrixValuedL2(mesh2d):
    mesh = mesh2d
    space = MatrixValued(L2(mesh, order=2), dim=2, symmetric=True)

    def _expr(x, y):
        return ((6 * x**2 + 4 * y**2, y**1 - 1.5 * x**2),
                (y**1 - 1.5 * x**2, y**2 - 1.5 * x**1))

    def _check(gf):
        return np.allclose(gf(mesh(0.5, 0.5)), np.array(_expr(0.5, 0.5)).flatten(), atol=1e-10, rtol=1e-10)

    _test_interpolation(space, _expr(x, y), _check)


def test_SymDevMatrixValuedL2(mesh2d):
    mesh = mesh2d
    space = MatrixValued(L2(mesh, order=2), dim=2, symmetric=True, deviatoric=True)

    def _expr(x, y):
        return ((6 * x**2 + 4 * y**2, y**1 - 1.5 * x**2),
                (y**1 - 1.5 * x**2, -6 * x**2 - 4 * y**2))

    def _check(gf):
        return np.allclose(gf(mesh(0.5, 0.5)), np.array(_expr(0.5, 0.5)).flatten(), atol=1e-10, rtol=1e-10)

    _test_interpolation(space, _expr(x, y), _check)


if __name__ == "__main__":
    mesh = mk_2d_mesh()
    test_MatrixValuedL2(mesh)
