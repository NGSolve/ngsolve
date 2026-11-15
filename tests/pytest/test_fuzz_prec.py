"""
Preconditioners: illegal use is rejected with a clean exception,
legal use (including creators without a bilinear form) keeps working.
"""
import pytest
from ngsolve import *
from ngsolve.meshes import MakeStructured2DMesh, MakeStructured3DMesh, Make1DMesh
from ngsolve.krylovspace import CGSolver
import ngsolve.comp


def poisson(mesh, order=1, **bfflags):
    fes = H1(mesh, order=order, dirichlet=".*")
    u, v = fes.TnT()
    a = BilinearForm(grad(u)*grad(v)*dx, **bfflags)
    f = LinearForm(v*dx)
    return fes, a, f


def cg_solve(a, f, pre, fes):
    a.Assemble()
    f.Assemble()
    inv = CGSolver(a.mat, pre, tol=1e-10, maxiter=200)
    gfu = GridFunction(fes)
    gfu.vec.data = inv * f.vec
    ref = GridFunction(fes)
    ref.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec
    diff = gfu.vec.CreateVector()
    diff.data = gfu.vec - ref.vec
    return inv.iterations, diff.Norm() / ref.vec.Norm()


@pytest.mark.parametrize("quads", [True, False])
@pytest.mark.parametrize("symmetric", [True, False])
def test_multigrid_skipped_levels(quads, symmetric):
    # assemble on level 0, refine twice, assemble again; refine before first assembly
    for assemble_first in [True, False]:
        m = MakeStructured2DMesh(nx=4, ny=4, quads=quads)
        fes, a, f = poisson(m, symmetric=symmetric, symmetric_storage=symmetric)
        c = Preconditioner(a, "multigrid")
        if assemble_first:
            a.Assemble()
        m.Refine(); m.Refine()
        its, err = cg_solve(a, f, c, fes)
        info = f"assemble_first={assemble_first}, inverse={a.mat.GetInverseType()}"
        assert err < 1e-8, info
        assert its < 25, info

    # the prolongation matrix agrees with the prolongation operator
    P = fes.Prolongation()
    for level in [1, 2]:
        Pmat = P.CreateMatrix(level)
        nc, nf = Pmat.width, Pmat.height
        v = BaseVector(nf); v[:] = 0
        v[0:nc] = 1
        P.Prolongate(level, v)
        w = BaseVector(nf)
        e = BaseVector(nc); e[:] = 1
        w.data = Pmat * e
        for i in range(nf):
            assert w[i] == pytest.approx(v[i])


@pytest.mark.parametrize("quads", [True, False])
@pytest.mark.parametrize("symmetric", [True, False])
def test_galerkin_projection_flags(quads, symmetric):
    # Galerkin coarse matrices equal the assembled ones and keep the symmetry flag
    # (pardiso reads symmetric storage from it)
    m = MakeStructured2DMesh(nx=4, ny=4, quads=quads)
    fes = FESpace("nodal", m, order=1, dirichlet=".*")
    u, v = fes.TnT()
    a = BilinearForm(grad(u)*grad(v)*dx, symmetric=symmetric, symmetric_storage=symmetric)
    c = Preconditioner(a, "multigrid")
    m.Refine(); m.Refine()
    a.Assemble()
    fine = a.GetMatrixLevel(2)
    for l in range(2):
        ml = MakeStructured2DMesh(nx=4, ny=4, quads=quads)
        for r in range(l):
            ml.Refine()
        fl = FESpace("nodal", ml, order=1)
        ul, vl = fl.TnT()
        ref = BilinearForm(grad(ul)*grad(vl)*dx, symmetric=symmetric, symmetric_storage=symmetric).Assemble().mat
        mat = a.GetMatrixLevel(l)
        assert type(mat) == type(ref)
        assert str(mat.is_symmetric) == str(fine.is_symmetric) == str(ref.is_symmetric)
        x = ref.CreateColVector()
        x.SetRandom()
        y = x.CreateVector()
        y.data = mat * x - ref * x
        assert y.Norm() < 1e-10 * x.Norm()


def test_precond_not_ready():
    m = MakeStructured2DMesh(nx=2, ny=2)
    with pytest.raises(Exception, match="not ready"):
        Preconditioner(BilinearForm(H1(m)), "local").ToDense()
    with pytest.raises(Exception, match="no bilinear form"):
        ngsolve.comp.H1AMG().height
    a = BilinearForm(H1(m), symmetric=True)
    c = ngsolve.comp.H1AMG(a)
    with pytest.raises(Exception, match="not ready"):
        c.height
    m3 = MakeStructured3DMesh(hexes=False, nx=1, ny=1, nz=1)
    for name in ["hcurlamg", "APhiamg"]:
        a3 = BilinearForm(HCurl(m3, order=0))
        c = Preconditioner(a3, name)
        with pytest.raises(Exception, match="not ready"):
            c.height
        with pytest.raises(Exception, match="not ready"):
            c.GetOperatorInfo()
    c = Preconditioner(BilinearForm(H1(m)), "multigrid")
    with pytest.raises(Exception, match="not ready"):
        c.height
    with pytest.raises(Exception, match="no bilinear form"):
        ngsolve.comp.BDDCPreconditioner().height


def test_bddc_created_after_assemble():
    m = MakeStructured2DMesh(nx=4, ny=4)
    fes, a, f = poisson(m, order=3)
    a.Assemble()
    c = Preconditioner(a, "bddc")
    w = a.mat.CreateColVector(); w[:] = 1
    y = w.CreateVector()
    with pytest.raises(Exception, match="before the bilinear form is assembled"):
        y.data = c * w
    with pytest.raises(Exception, match="before assembling"):
        c.Update()
    # a fresh assembly sets it up
    its, err = cg_solve(a.Assemble(reallocate=True), f, c, fes)
    assert err < 1e-8 and its < 20


def test_bddc_scalar_type():
    m = MakeStructured2DMesh(nx=4, ny=4)
    fes, a, f = poisson(m, order=2)
    for name in ["bddcc", "bddcrc"]:
        with pytest.raises(Exception, match="does not fit"):
            Preconditioner(a, name)
    # bddc picks the type from the form, also for real matrices on a complex space
    fesc = H1(m, order=3, dirichlet=".*", complex=True)
    u, v = fesc.TnT()
    for flags in [{}, {"real": True}]:
        ac = BilinearForm(grad(u)*grad(v)*dx, **flags)
        fc = LinearForm(1j*v*dx)
        c = Preconditioner(ac, "bddc")
        its, err = cg_solve(ac, fc, c, fesc)
        assert err < 1e-8 and its < 20
    # wrong vectors for bddcrc are rejected by the operator check
    ac = BilinearForm(grad(u)*grad(v)*dx, real=True)
    c = Preconditioner(ac, "bddcrc")
    ac.Assemble()
    assert c.is_complex
    x = BaseVector(fesc.ndof, complex=False); x[:] = 1
    y = x.CreateVector()
    with pytest.raises(Exception):
        y.data = c * x


def test_multigrid_creator():
    m = MakeStructured2DMesh(nx=4, ny=4, quads=False)
    fes, a, f = poisson(m)
    coarse = BilinearForm(fes)
    creator = ngsolve.comp.MultiGridPreconditioner(lo_preconditioner=IdentityMatrix(fes.ndof))
    assert creator.IsCreator()
    with pytest.raises(Exception, match="no bilinear form"):
        creator.height
    with pytest.raises(Exception, match="no bilinear form"):
        creator.Update()
    str(creator)
    with pytest.raises(Exception, match="no bilinear form"):
        BlockMatrix([[ngsolve.comp.MultiGridPreconditioner(name="x")]]).ToDense()

    # creators are legal and create working preconditioners
    for creator in [ngsolve.comp.MultiGridPreconditioner(), ngsolve.comp.H1AMG(),
                    ngsolve.comp.BDDCPreconditioner(), ngsolve.comp.LocalPreconditioner()]:
        fes, a, f = poisson(m)
        c = creator.Create(a)
        its, err = cg_solve(a, f, c, fes)
        assert err < 1e-8 and its < 60


def test_amg_wrong_space():
    m2q = MakeStructured2DMesh(nx=2, ny=2)
    m2t = MakeStructured2DMesh(nx=2, ny=2, quads=False)
    m3 = MakeStructured3DMesh(hexes=False, nx=2, ny=2, nz=2)

    def assemble(fes, form, name):
        u, v = fes.TnT()
        a = BilinearForm(form(u, v))
        c = Preconditioner(a, name)
        a.Assemble()
        return a, c

    curlcurl = lambda u, v: curl(u)*curl(v)*dx + u*v*dx
    with pytest.raises(Exception, match="triangular faces"):
        assemble(HCurl(m2q, order=0, complex=True), curlcurl, "hcurlamg")
    with pytest.raises(Exception, match="lowest order HCurl"):
        assemble(HCurlCurl(m2t, order=0), lambda u, v: InnerProduct(u, v)*dx, "hcurlamg")
    with pytest.raises(Exception, match="APhiamg needs"):
        assemble(L2(Make1DMesh(4), order=9), lambda u, v: u*v*dx, "APhiamg")
    with pytest.raises(Exception, match="APhiamg needs"):
        assemble(HDiv(m3, order=0), lambda u, v: div(u)*div(v)*dx + u*v*dx, "APhiamg")
    with pytest.raises(Exception, match="does not fit"):
        ngsolve.comp.H1AMG(BilinearForm(H1(m2t, complex=True)))

    # legal: real and complex lowest order HCurl on simplices, complex h1amg
    for fes in [HCurl(m2t, order=0, dirichlet=".*"), HCurl(m3, order=0, complex=True, dirichlet=".*")]:
        u, v = fes.TnT()
        a = BilinearForm(curl(u)*curl(v)*dx + u*v*dx)
        f = LinearForm(fes)
        f += (CF((1, 0)) if fes.mesh.dim == 2 else CF((1, 0, 0))) * v * dx
        c = Preconditioner(a, "hcurlamg")
        its, err = cg_solve(a, f, c, fes)
        assert err < 1e-8 and its < 40
    fes = H1(m3, order=1, complex=True, dirichlet=".*")
    u, v = fes.TnT()
    a = BilinearForm(grad(u)*grad(v)*dx)
    f = LinearForm(1j*v*dx)
    c = Preconditioner(a, "h1amg")
    its, err = cg_solve(a, f, c, fes)
    assert err < 1e-8 and its < 30
