import pytest
from ngsolve import *
from ngsolve.solvers import mkl_pardiso

pytestmark = pytest.mark.skipif(mkl_pardiso._pardiso is None, reason="MKL Pardiso not available")


@pytest.mark.parametrize("is_complex", [False, True])
@pytest.mark.parametrize("symmetric", [False, True])
def test_pardiso(is_complex, symmetric):
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.1))
    fes = H1(mesh, order=2, dirichlet="left", complex=is_complex)
    u, v = fes.TnT()
    form = grad(u) * grad(v) * dx + (1j if is_complex else 1) * u * v * dx
    if not symmetric:
        form += 0.3 * grad(u)[0] * v * dx
    a = BilinearForm(form, symmetric=symmetric).Assemble()
    f = LinearForm(x * v * dx).Assemble()

    inv = a.mat.Inverse(fes.FreeDofs(), inverse="pardiso")
    proj = Projector(fes.FreeDofs(), True)
    u = f.vec.CreateVector()
    for op in [inv, inv.T]:
        mat = a.mat if op is inv else a.mat.T
        u.data = op * f.vec
        assert Norm(proj * (f.vec - mat * u)) < 1e-10 * Norm(f.vec)


def test_pardiso_temporary_matrix():
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.2))
    fes = H1(mesh, order=1)
    u, v = fes.TnT()
    inv = BilinearForm(u * v * dx).Assemble().mat.Inverse(inverse="pardiso")
    f = LinearForm(v * dx).Assemble()
    x = f.vec.CreateVector()
    x.data = inv * f.vec
    assert Norm(x) > 0


def test_pardiso_empty_and_complex_vectors():
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.2))
    # all dofs local: nothing left after condensation
    fes = L2(mesh, order=1)
    u, v = fes.TnT()
    a = BilinearForm(u * v * dx, condense=True).Assemble()
    x = a.mat.CreateVector()
    x[:] = 1
    y = x.CreateVector()
    y.data = a.mat.Inverse(fes.FreeDofs(True), inverse="pardiso") * x
    assert Norm(y) == 0

    # real matrix on a complex space
    fes = H1(mesh, order=2, dirichlet="left", complex=True)
    u, v = fes.TnT()
    a = BilinearForm(grad(u) * grad(v) * dx, real=True).Assemble()
    f = LinearForm(1j * v * dx).Assemble()
    inv = a.mat.Inverse(fes.FreeDofs(), inverse="pardiso")
    gfu = GridFunction(fes)
    gfu.vec.data = inv * f.vec
    proj = Projector(fes.FreeDofs(), True)
    assert Norm(proj * (f.vec - a.mat * gfu.vec)) < 1e-10 * Norm(f.vec)
