import pytest
from ngsolve import *
from ngsolve.comp import RegisterPreconditioner
import netgen.occ as occ

# Dirichlet conditions of the compact Solve must reach the preconditioner

@pytest.fixture(scope="module")
def mesh():
    box = occ.Box((0,0,0), (1,1,1))
    box.faces.Min(occ.Z).name = "bot"
    box.faces.Max(occ.Z).name = "top"
    return Mesh(occ.OCCGeometry(box).GenerateMesh(maxh=0.3))

def cg():
    return solvers.CGSolver.Creator(maxiter=2000, tol=1e-12)

def err(gf, ref):
    return sqrt(abs(Integrate(InnerProduct(gf-ref, gf-ref), gf.space.mesh)))

def solve_with(mkpre, lhs, rhs, bc):
    # as Solve does, for preconditioners without creator mode
    a = BilinearForm(lhs-rhs)
    pre = mkpre(a)
    gf = GridFunction(a.space)
    gf[bc.proxy][bc.vbn] = bc.val
    a.AssembleLinearization(gf.vec)
    inv = solvers.CGSolver(a.mat, pre, maxiter=2000, tol=1e-12)
    gf.vec.data -= inv * a.Apply(gf.vec)
    return gf


@pytest.mark.parametrize("is_complex", [False, True])
@pytest.mark.parametrize("val, expected", [(0, 0), (1+x, 1.5)])
def test_h1amg(mesh, is_complex, val, expected):
    fes = H1(mesh, order=1, complex=is_complex)
    u, v = fes.TnT()
    lhs, rhs = grad(u)*grad(v)*dx + u*v*ds("top"), 1*v*dx
    bc = u[BND("bot")] == val
    ref = Solve(lhs == rhs, bc)
    gf = Solve(lhs == rhs, bc, preconditioners.H1AMG(), cg())
    assert err(gf, ref) < 1e-10
    assert gf(mesh(0.5, 0.5, 0)) == pytest.approx(expected, abs=1e-10)
    gf = solve_with(lambda a: Preconditioner(a, "h1amg", additional_dirbc=[bc]), lhs, rhs, bc)
    assert err(gf, ref) < 1e-10


def test_hcurlamg(mesh):
    fes = HCurl(mesh, order=0)
    u, v = fes.TnT()
    lhs, rhs = curl(u)*curl(v)*dx + 0.1*u*v*dx, CF((1,y,0))*v*dx
    bc = u[BND("bot")] == 0
    ref = Solve(lhs == rhs, bc)
    gf = solve_with(lambda a: preconditioners.HCurlAMG(a, additional_dirbc=[bc]), lhs, rhs, bc)
    assert err(gf, ref) < 1e-9


def test_python_preconditioner(mesh):
    seen = []
    def smoother(mat, freedofs, *args, **kwargs):
        seen.append(freedofs.NumSet())
        return mat.CreateSmoother(freedofs)
    RegisterPreconditioner("test_dirbc_smoother", smoother)

    fes = H1(mesh, order=1)
    u, v = fes.TnT()
    lhs, rhs = grad(u)*grad(v)*dx + u*v*ds("top"), 1*v*dx
    bc = u[BND("bot")] == 0
    ref = Solve(lhs == rhs, bc)
    gf = solve_with(lambda a: Preconditioner(a, "test_dirbc_smoother", additional_dirbc=[bc]), lhs, rhs, bc)
    assert seen[-1] == fes.ndof - fes.GetDofs(mesh.Boundaries("bot")).NumSet()
    assert err(gf, ref) < 1e-10
