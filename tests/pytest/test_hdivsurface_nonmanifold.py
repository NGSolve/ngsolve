import math
import pytest
from ngsolve import *
from netgen.occ import *

# sheets glued along the z-axis: non-manifold junction for 3 and more sheets,
# for 2 sheets a manifold edge where the faces are inconsistently oriented

def SheetMesh(nsheets, maxh=0.3, angles=None):
    angles = angles or [2*math.pi*k/nsheets for k in range(nsheets)]
    faces = []
    for k, th in enumerate(angles):
        d = (math.cos(th), math.sin(th), 0)
        p = [Pnt(0,0,0), Pnt(d[0],d[1],0), Pnt(d[0],d[1],1), Pnt(0,0,1)]
        f = Face(Wire([Segment(p[i], p[(i+1)%4]) for i in range(4)]))
        f.name = f"sheet{k}"
        faces.append(f)
    shape = Glue(faces)
    for e in shape.edges:
        c = e.center
        e.name = "junction" if abs(c.x) < 1e-8 and abs(c.y) < 1e-8 else "outer"
    return Mesh(OCCGeometry(shape).GenerateMesh(maxh=maxh)), angles


def MixedPoisson(mesh, order, condense=False):
    # -Laplace u = 1, u = 0 on outer edges, compared to a high order primal solution
    Vh = H1(mesh, order=order+3, dirichlet_bbnd="outer")
    u, v = Vh.TnT()
    a = BilinearForm(grad(u).Trace()*grad(v).Trace()*ds).Assemble()
    l = LinearForm(v*ds).Assemble()
    uref = GridFunction(Vh)
    uref.vec.data = a.mat.Inverse(Vh.FreeDofs(), "sparsecholesky") * l.vec

    X = HDivSurface(mesh, order=order) * SurfaceL2(mesh, order=order-1)
    (sig, uu), (tau, vv) = X.TnT()
    am = BilinearForm(X, condense=condense)
    am += (sig.Trace()*tau.Trace() + uu*div(tau.Trace()) + div(sig.Trace())*vv)*ds
    am.Assemble()
    lm = LinearForm(-vv*ds).Assemble()
    gfx = GridFunction(X)
    if condense:
        fr = lm.vec.CreateVector()
        fr.data = lm.vec + am.harmonic_extension_trans * lm.vec
        gfx.vec.data = am.mat.Inverse(X.FreeDofs(True), "umfpack") * fr
        gfx.vec.data += am.inner_solve * lm.vec
        gfx.vec.data += am.harmonic_extension * gfx.vec
    else:
        gfx.vec.data = am.mat.Inverse(X.FreeDofs(), "umfpack") * lm.vec
    gsig, gu = gfx.components
    erru = sqrt(Integrate((gu-uref)**2, mesh, BND)) / sqrt(Integrate(uref**2, mesh, BND))
    errs = sqrt(Integrate((gsig.Trace()-grad(uref).Trace())**2, mesh, BND))
    return erru, errs


@pytest.mark.parametrize("nsheets", [2, 3, 4])
def test_mixed_poisson(nsheets):
    mesh, _ = SheetMesh(nsheets, angles=[0, math.pi/2] if nsheets == 2 else None)
    erru, errs = MixedPoisson(mesh, 2)
    # with one dof per edge mode at the junction the error was about 50 percent
    assert erru < 0.05
    assert errs < 0.01


def test_mixed_poisson_condense():
    mesh, _ = SheetMesh(3)
    assert MixedPoisson(mesh, 2) == pytest.approx(MixedPoisson(mesh, 2, condense=True), rel=1e-8)


@pytest.mark.parametrize("order", [1, 2])
def test_set_junction_flux(order):
    # constant flux on every sheet, the fluxes into the junction sum up to zero
    cs = [1, -2, 0.5, 0.5]
    mesh, angles = SheetMesh(4, maxh=0.4, angles=[0, 1.0, 2.5, 4.0])
    cf = mesh.BoundaryCF({f"sheet{k}": CF((cs[k]*math.cos(th), cs[k]*math.sin(th), 0))
                          for k, th in enumerate(angles)})
    gf = GridFunction(HDivSurface(mesh, order=order))
    gf.Set(cf, definedon=mesh.Boundaries(".*"))
    assert sqrt(Integrate((gf.Trace()-cf)**2, mesh, BND)) < 1e-12


def test_apply_energy_linearization():
    # matrix-free operations agree with the assembled matrix
    mesh, _ = SheetMesh(3, maxh=0.4)
    fes = HDivSurface(mesh, order=2)
    s, t = fes.TnT()
    x = GridFunction(fes).vec
    x.SetRandom()
    a = BilinearForm((s.Trace()*t.Trace() + div(s.Trace())*div(t.Trace()))*ds).Assemble()
    y = x.CreateVector()
    a.Apply(x, y)
    assert Norm(y - a.mat*x) < 1e-10 * Norm(y)
    assert a.Energy(x) == pytest.approx(0.5*InnerProduct(x, a.mat*x), rel=1e-12)

    an = BilinearForm((s.Trace()*s.Trace())*(s.Trace()*t.Trace())*ds)
    an.AssembleLinearization(x)
    d = x.CreateVector()
    d.SetRandom()
    eps = 1e-6
    x2 = x.CreateVector()
    x2.data = x + eps*d
    y1, y2 = x.CreateVector(), x.CreateVector()
    an.Apply(x, y1)
    an.Apply(x2, y2)
    fd = (1/eps) * (y2 - y1).Evaluate()
    assert Norm(fd - an.mat*d) < 1e-4 * Norm(fd)
