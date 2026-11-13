from ngsolve import *
from netgen import occ
import pytest

# two rectangles with a non-matching interface at y = 0.5
@pytest.fixture(scope="module")
def mesh():
    r1 = occ.MoveTo(0,0).Rectangle(1,0.5).Face()
    r2 = occ.MoveTo(0,0.5).Rectangle(1,0.5).Face()
    r1.edges.Max(occ.Y).name = "top1"
    r2.edges.Min(occ.Y).name = "bot2"
    r1.edges.Min(occ.Y).name = "fix"
    r2.edges.Max(occ.Y).name = "fix"
    r1.faces.name = "lower"
    r2.faces.name = "upper"
    r2.faces.maxh = 0.07
    geo = occ.OCCGeometry(occ.Compound([r1,r2]), dim=2)
    return Mesh(geo.GenerateMesh(maxh=0.1))

def Penalty(u, v):
    return 20/specialcf.mesh_size*(u-u.Other())*(v-v.Other())

def Nitsche(u, v):
    n = specialcf.normal(2)
    return Penalty(u,v) + n*grad(u)*(v.Other()-v) + n*grad(v)*(u.Other()-u)

def TestVector(fes):
    gf = GridFunction(fes)
    gf.Set(x*y+sin(3*x)+IfPos(y-0.5, x*x, 0))
    return gf.vec

def Shift(mesh, dx):
    defo = GridFunction(VectorH1(mesh, order=1))
    defo.Set(CF((dx,0)), definedon=mesh.Materials("upper"))
    return defo

def Diff(a, b, vec):
    return Norm((a.mat*vec - b.mat*vec).Evaluate())

# reference: the same integrals by ContactBoundary
def Reference(mesh, a, form, defo=None, volume=False, **kwargs):
    cb = ContactBoundary(mesh.Boundaries("top1"), mesh.Boundaries("bot2"), volume=volume)
    cb.AddIntegrator(form)
    cb.Update(defo, bf=a, **kwargs)
    return cb


@pytest.mark.parametrize("quadrature", ["gauss", "clip"])
@pytest.mark.parametrize("regions", [False, True])
def test_assemble(mesh, quadrature, regions):
    fes = H1(mesh, order=2, dirichlet="fix")
    u,v = fes.TnT()
    me, other = "top1", "bot2"
    if regions:
        me, other = mesh.Boundaries(me), mesh.Boundaries(other)
    a = BilinearForm(grad(u)*grad(v)*dx)
    a += ContactIntegrator(me, other, intorder=4, quadrature=quadrature).Add(Penalty(u,v))
    a.Assemble()
    b = BilinearForm(grad(u)*grad(v)*dx)
    cb = Reference(mesh, b, Penalty(u,v), intorder=4, quadrature=quadrature)
    b.Assemble()
    plain = BilinearForm(grad(u)*grad(v)*dx).Assemble()
    vec = TestVector(fes)
    assert Diff(a, plain, vec) > 1
    assert Diff(a, b, vec) < 1e-12


def test_volume_terms(mesh):
    fes = H1(mesh, order=2, dirichlet="fix")
    u,v = fes.TnT()
    a = BilinearForm(grad(u)*grad(v)*dx)
    a += ContactIntegrator("top1", "bot2", intorder=6, volume=True).Add(Nitsche(u,v))
    a.Assemble()
    b = BilinearForm(grad(u)*grad(v)*dx)
    cb = Reference(mesh, b, Nitsche(u,v), volume=True, intorder=6)
    b.Assemble()
    assert Diff(a, b, TestVector(fes)) < 1e-12


# Nitsche coupling with exact quadrature reproduces a linear solution
def test_nitsche_consistency(mesh):
    fes = H1(mesh, order=2, dirichlet="fix")
    u,v = fes.TnT()
    a = BilinearForm(grad(u)*grad(v)*dx)
    a += ContactIntegrator("top1", "bot2", intorder=4, quadrature="clip", volume=True).Add(Nitsche(u,v))
    a.Assemble()
    gfu = GridFunction(fes)
    gfu.Set(y, BND)
    res = (-a.mat * gfu.vec).Evaluate()
    gfu.vec.data += a.mat.Inverse(fes.FreeDofs()) * res
    assert Integrate((gfu-y)**2*dx, mesh) < 1e-20


def test_reassemble(mesh):
    fes = H1(mesh, order=2, dirichlet="fix")
    u,v = fes.TnT()
    defo = Shift(mesh, 0)
    a = BilinearForm(grad(u)*grad(v)*dx)
    a += ContactIntegrator("top1", "bot2", deformation=defo, intorder=4, maxdist=0.1).Add(Penalty(u,v))
    a.Assemble()
    nze = a.mat.nze
    # Assemble searches the pairs again
    defo.vec.data = Shift(mesh, 0.33).vec
    a.Assemble()
    assert a.mat.nze != nze
    b = BilinearForm(grad(u)*grad(v)*dx)
    cb = Reference(mesh, b, Penalty(u,v), defo, intorder=4, maxdist=0.1)
    b.Assemble()
    assert Diff(a, b, TestVector(fes)) < 1e-12


def test_apply_and_linearization(mesh):
    fes = H1(mesh, order=2, dirichlet="fix")
    u,v = fes.TnT()
    a = BilinearForm(grad(u)*grad(v)*dx)
    a += ContactIntegrator("top1", "bot2", intorder=4).Add(Penalty(u,v))
    a.Assemble()
    vec = TestVector(fes)
    ax = (a.mat * vec).Evaluate()
    assert Norm((a.Apply(vec) - ax).Evaluate()) < 1e-10 * Norm(ax)
    a.AssembleLinearization(vec)
    assert Norm((a.mat * vec - ax).Evaluate()) < 1e-10 * Norm(ax)


def test_apply_without_assemble(mesh):
    fes = H1(mesh, order=2, dirichlet="fix")
    u,v = fes.TnT()
    defo = Shift(mesh, 0)
    a = BilinearForm(fes)
    ci = ContactIntegrator("top1", "bot2", deformation=defo, intorder=4, maxdist=0.1).Add(Penalty(u,v))
    a += ci
    vec = TestVector(fes)
    y0 = a.Apply(vec).Evaluate()
    assert Norm(y0) > 1
    # the pairs are kept until Update is called
    defo.vec.data = Shift(mesh, 0.1).vec
    assert Norm((a.Apply(vec) - y0).Evaluate()) == 0
    ci.Update()
    y1 = a.Apply(vec).Evaluate()
    assert Norm((y1 - y0).Evaluate()) > 1e-3
    b = BilinearForm(fes)
    cb = Reference(mesh, b, Penalty(u,v), defo, intorder=4, maxdist=0.1)
    assert Norm((b.Apply(vec) - y1).Evaluate()) < 1e-12


def test_mixed_spaces(mesh):
    fes = H1(mesh, order=2, dirichlet="fix")
    fes2 = H1(mesh, order=1, dirichlet="fix")
    u,v = fes.TrialFunction(), fes2.TestFunction()
    a = BilinearForm(trialspace=fes, testspace=fes2)
    a += grad(u)*grad(v)*dx
    a += ContactIntegrator("top1", "bot2", intorder=4).Add(Penalty(u,v))
    a.Assemble()
    b = BilinearForm(trialspace=fes, testspace=fes2)
    b += grad(u)*grad(v)*dx
    cb = Reference(mesh, b, Penalty(u,v), intorder=4)
    b.Assemble()
    vec = TestVector(fes)
    assert a.mat.nze == b.mat.nze
    assert Diff(a, b, vec) < 1e-12
    assert Norm((a.Apply(vec) - a.mat*vec).Evaluate()) < 1e-10


def test_energy(mesh):
    fes = VectorH1(mesh, order=2, dirichlet="fix")
    u = fes.TrialFunction()
    gfu = GridFunction(fes)
    ci = ContactIntegrator(mesh.Boundaries("top1"), mesh.Boundaries("bot2"),
                           deformation=gfu, intorder=4, maxdist=0.1)
    X = CF((x,y))
    cf = (X + u - X.Other() - u.Other()) * ci.normal
    ci.AddEnergy(IfPos(cf, 1e6*cf*cf, 0))
    a = BilinearForm(fes)
    a += ci

    def Move(gf, val):
        gf.Set(CF((0,val)), definedon=mesh.Materials("upper"))
    # separated
    Move(gfu, 0.01)
    assert a.Energy(gfu.vec) == 0
    # penetrating
    Move(gfu, -0.01)
    energy = a.Energy(gfu.vec)
    assert energy == pytest.approx(1e6*0.01**2, rel=1e-8)

    # Apply is the derivative of the energy, the linearization the derivative of Apply
    w = GridFunction(fes)
    Move(w, 1)
    eps = 1e-5
    ur = (gfu.vec + eps*w.vec).Evaluate()
    ul = (gfu.vec - eps*w.vec).Evaluate()
    d = a.Apply(gfu.vec).Evaluate()
    assert (a.Energy(ur) - a.Energy(ul))/2/eps == pytest.approx(InnerProduct(d, w.vec), rel=1e-8)
    a.AssembleLinearization(gfu.vec)
    dd = (1/(2*eps) * (a.Apply(ur) - a.Apply(ul))).Evaluate()
    assert Norm((a.mat*w.vec - dd).Evaluate()) < 1e-6 * Norm(dd)


def test_checks(mesh):
    fes = H1(mesh, order=2, dirichlet="fix")
    fes1 = H1(mesh, order=1)
    u,v = fes.TnT()
    u1,v1 = fes1.TnT()

    # normal needs the mesh
    ci = ContactIntegrator("top1", "bot2")
    with pytest.raises(Exception):
        ci.normal
    ci.Add(Penalty(u,v))
    assert ci.normal.dim == 2
    assert ContactIntegrator(mesh.Boundaries("top1"), mesh.Boundaries("bot2")).normal.dim == 2

    with pytest.raises(Exception):
        ci.Add(Penalty(u1,v1))

    # form on another space than the BilinearForm
    a = BilinearForm(grad(u)*grad(v)*dx)
    a += ContactIntegrator("top1", "bot2").Add(Penalty(u1,v1))
    with pytest.raises(Exception):
        a.Assemble()

    # nothing added
    a = BilinearForm(grad(u)*grad(v)*dx)
    a += ContactIntegrator("top1", "bot2")
    a.Assemble()
    plain = BilinearForm(grad(u)*grad(v)*dx).Assemble()
    assert Diff(a, plain, TestVector(fes)) == 0
