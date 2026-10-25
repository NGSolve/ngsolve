from ngsolve import *
from netgen.geom2d import unit_square
from netgen.csg import unit_cube
import numpy as np
import pytest


def MakeMesh(dim):
    if dim == 2:
        return Mesh(unit_square.GenerateMesh(maxh=0.3))
    return Mesh(unit_cube.GenerateMesh(maxh=0.4))


def RelDiff(a, b):
    return Norm(a.mat.AsVector()-b.mat.AsVector()) / Norm(a.mat.AsVector())


@pytest.mark.parametrize("dim", [2, 3])
def test_strain_values(dim):
    mesh = MakeMesh(dim)
    fes = VectorH1(mesh, order=3)
    gf = GridFunction(fes)
    gf.Set(CF(tuple(sin(x+2*(i+1)*y) + x*y*(i+1) for i in range(dim))))
    eps = gf.Operator("strain")
    assert tuple(eps.dims) == (dim, dim)
    diff = eps - Sym(Grad(gf))
    assert sqrt(Integrate(InnerProduct(diff, diff), mesh)) < 1e-12


@pytest.mark.parametrize("dim", [2, 3])
@pytest.mark.parametrize("interleaved", [False, True])
def test_strain_matrices(dim, interleaved):
    mesh = MakeMesh(dim)
    fes = VectorH1(mesh, order=2, interleaved=interleaved)
    u, v = fes.TnT()
    eu, ev = u.Operator("strain"), v.Operator("strain")
    C = CF(tuple(range(1, dim*dim+1)), dims=(dim, dim))
    forms = [
        (InnerProduct(Sym(Grad(u)), Sym(Grad(v))), InnerProduct(eu, ev)),
        (Trace(Sym(Grad(u)))*Trace(Sym(Grad(v))), Trace(eu)*Trace(ev)),
        (2*InnerProduct(Sym(Grad(u)), Sym(Grad(v))) + div(u)*div(v),
         2*InnerProduct(eu, ev) + Trace(eu)*Trace(ev)),
        (InnerProduct(C*Sym(Grad(u)), Sym(Grad(v))), InnerProduct(C*eu, ev)),
        (InnerProduct(Grad(u), Sym(Grad(v))), InnerProduct(Grad(u), ev)),
    ]
    for f_ref, f_strain in forms:
        a_ref = BilinearForm(f_ref*dx).Assemble()
        a_strain = BilinearForm(f_strain*dx).Assemble()
        assert RelDiff(a_ref, a_strain) < 1e-12


@pytest.mark.parametrize("dim", [2, 3])
def test_strain_linearform_and_apply(dim):
    mesh = MakeMesh(dim)
    fes = VectorH1(mesh, order=2)
    u, v = fes.TnT()
    sigma = CF(tuple(x*(i+1)+y for i in range(dim*dim)), dims=(dim, dim))
    f_ref = LinearForm(InnerProduct(sigma, Sym(Grad(v)))*dx).Assemble()
    f_strain = LinearForm(InnerProduct(sigma, v.Operator("strain"))*dx).Assemble()
    assert Norm(f_ref.vec-f_strain.vec) / Norm(f_ref.vec) < 1e-12

    eu, ev = u.Operator("strain"), v.Operator("strain")
    a = BilinearForm(InnerProduct(eu, ev)*dx).Assemble()
    afree = BilinearForm(InnerProduct(eu, ev)*dx, nonassemble=True)
    vx = a.mat.CreateColVector()
    vx.FV().NumPy()[:] = np.sin(np.arange(len(vx)))
    vy = vx.CreateVector()
    afree.Apply(vx, vy)
    vy -= a.mat * vx
    assert Norm(vy) / Norm(vx) < 1e-12


# matrix-valued spaces with vector-space embedding: x^T A y must equal the integral
@pytest.mark.parametrize("kwargs", [dict(symmetric=True), dict(skewsymmetric=True),
                                    dict(symmetric=True, deviatoric=True)])
def test_matrixvalued_embedding(kwargs):
    mesh = MakeMesh(3)
    fes = MatrixValued(H1(mesh, order=2), **kwargs)
    u, v = fes.TnT()
    gx, gy = GridFunction(fes), GridFunction(fes)
    gx.vec.FV().NumPy()[:] = np.sin(np.arange(fes.ndof))
    gy.vec.FV().NumPy()[:] = np.cos(np.arange(fes.ndof))
    C = CF(tuple(0.1*np.arange(1, 10) + np.eye(3).flatten()), dims=(3, 3))
    forms = [
        (InnerProduct(u, v), InnerProduct(gx, gy)),
        (InnerProduct(C*u, v), InnerProduct(C*gx, gy)),
        (Trace(u)*Trace(v) + InnerProduct(u, v), Trace(gx)*Trace(gy) + InnerProduct(gx, gy)),
    ]
    for f, fg in forms:
        a = BilinearForm(f*dx).Assemble()
        lhs = InnerProduct(gy.vec, a.mat*gx.vec)
        rhs = Integrate(fg, mesh, order=8)
        assert abs(lhs-rhs) <= 1e-10 * abs(rhs)


def test_strain_complex_coefficient():
    mesh = MakeMesh(3)
    fes = VectorH1(mesh, order=2, complex=True)
    u, v = fes.TnT()
    eu, ev = u.Operator("strain"), v.Operator("strain")
    a_ref = BilinearForm((2+1j)*InnerProduct(Sym(Grad(u)), Sym(Grad(v)))*dx
                         + 1j*div(u)*div(v)*dx).Assemble()
    a_strain = BilinearForm((2+1j)*InnerProduct(eu, ev)*dx
                            + 1j*Trace(eu)*Trace(ev)*dx).Assemble()
    assert RelDiff(a_ref, a_strain) < 1e-12
