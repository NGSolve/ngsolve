from ngsolve import *
from netgen.csg import *
from netgen.geom2d import SplineGeometry
import pytest

@pytest.fixture
def mesh():
    geo = SplineGeometry()
    geo.AddRectangle( (-10, 0.0), (10, 1), leftdomain=1, bcs=["bottom","right","minion","left"])

    geo.AddCircle ( (-1, 5), r=1, leftdomain=2, bc="master")
    geo.SetMaterial(1, "brick")
    geo.SetMaterial(2, "ball")
    mesh = Mesh( geo.GenerateMesh(maxh=0.5))
    mesh.Curve(5)
    Draw(mesh)
    return mesh

def GetForms(fes):
    mesh = fes.mesh
    cb = ContactBoundary(mesh.Boundaries("master"), mesh.Boundaries("minion"))
    u = fes.TrialFunction()
    X = CoefficientFunction((x,y))
    cf = (u - u.Other() + X - X.Other()) * cb.normal
    cb.AddEnergy(IfPos(cf, 1e6 * cf * cf, 0))
    a = BilinearForm(fes)
    u = GridFunction(fes)
    return cb, a, u

def SetY(u, val):
    u.Set((0, val), definedon=u.space.mesh.Materials("ball"))

tested_spaces = [(VectorH1,{"order" : 3}), (H1, {"dim" : 2,
                                                 "order" : 3})]

@pytest.mark.parametrize("space, space_args", tested_spaces)
def test_energy(mesh, space, space_args):
    fes = space(mesh, **space_args)
    cb, a, u = GetForms(fes)
    SetY(u, -2.9)
    cb.Update(u, a, 4, 2)
    assert a.Energy(u.vec) == 0.
    SetY(u, -3.1)
    assert a.Energy(u.vec) > 1000
    SetY(u, -1)
    cb.Update(u, a, 4, 1)
    assert a.Energy(u.vec) == 0.
    SetY(u, -3.1)
    assert a.Energy(u.vec) == 0.

@pytest.mark.parametrize("space, space_args", tested_spaces)
def test_apply(mesh, space, space_args):
    fes = space(mesh, **space_args)
    cb, a, u = GetForms(fes)
    SetY(u, -2.9)
    cb.Update(u, a, 4, 2)
    SetY(u, -3.1)
    w = GridFunction(fes)
    SetY(w, -1)
    ur = u.vec.CreateVector()
    ul = u.vec.CreateVector()
    eps = 1e-4
    ur.data = u.vec + eps * w.vec
    ul.data = u.vec - eps * w.vec
    d = u.vec.CreateVector()
    a.Apply(u.vec, d)
    assert (a.Energy(ur) - a.Energy(ul))/2/eps == pytest.approx(InnerProduct(d,w.vec), rel=1e-10)

@pytest.mark.parametrize("space, space_args", tested_spaces)
def test_calclinearized(mesh, space, space_args):
    fes = space(mesh, **space_args)
    cb, a, u = GetForms(fes)
    SetY(u, -2.9)
    cb.Update(u, a, 4, 2)
    SetY(u, -3.1)
    w = GridFunction(fes)
    SetY(w, -1)
    ur = u.vec.CreateVector()
    ul = u.vec.CreateVector()
    eps = 1e-4
    ur.data = u.vec + eps * w.vec
    ul.data = u.vec - eps * w.vec
    d = u.vec.CreateVector()
    a.Apply(u.vec, d)
    vl = d.CreateVector()
    vr = d.CreateVector()
    a.Apply(ul, vl)
    a.Apply(ur, vr)
    a.AssembleLinearization(u.vec)
    d.data = a.mat * w.vec
    assert (InnerProduct(vr, w.vec)-InnerProduct(vl, w.vec))/2/eps == pytest.approx(InnerProduct(d, w.vec), rel=1e-10)

def test_gapfunction():
    geo = CSGeometry()
    r = 0.01
    center = (0.03, 0.02, 0.021)
    brick = OrthoBrick(Pnt(0,0,0), Pnt(0.1,0.04, 0.01)).bc("brick").mat("brick")
    ball = Sphere(Pnt(*center), r).bc("ball").mat("ball")
    geo.Add(brick)
    geo.Add(ball)
    mesh = Mesh( geo.GenerateMesh(maxh=0.01))

    mesh.Curve(order=3)
    Draw (mesh)

    master = Region(mesh, BND, "brick")
    minion = Region(mesh, BND, "ball")

    fes = H1(mesh, dim=mesh.dim)
    gfu = GridFunction(fes)
    gfu.vec[:] = 0.0

    cb = ContactBoundary(master, minion)
    cb.Update(gfu, maxdist=1.)
    Draw(cb.gap, mesh, "gap")
    Draw(Norm(cb.gap), mesh, "dist")

    error = Norm(-cb.gap + center - (x,y,z)) - r
    assert Integrate(error, mesh, definedon=master, order=1) < 1e-8

