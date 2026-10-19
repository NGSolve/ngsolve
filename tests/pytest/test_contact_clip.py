import numpy as np
import pytest
from ngsolve import *
from ngsolve.comp import ContactBoundary
from netgen import occ


@pytest.fixture(scope="module")
def split_cylinder():
    cyl = occ.Cylinder(occ.Pnt(0, 0, 0), occ.Vec(1, 0, 0), 1.0, 2.0)
    cut = occ.HalfSpace((1, 0, 0), (1, 0, 0))
    cut.faces.name = "interface"
    cyl.faces.name = "wall"
    right, left = cyl - cut, cyl * cut
    right.solids.name = "stat"
    left.solids.name = "rot"
    for part, tag in ((right, "_b"), (left, "_a")):
        for f in part.faces:
            f.name = f.name + tag
    geo = occ.Compound([right, left])
    geo.faces["interface_a"][0].Identify(geo.faces["interface_b"][0], "interface",
                                         occ.IdentificationType.PERIODIC)
    return Mesh(occ.OCCGeometry(geo).GenerateMesh(maxh=0.4))


@pytest.fixture(scope="module")
def split_rectangle():
    # two squares touching at x=1, meshed independently: non-matching interface
    left = occ.WorkPlane().Rectangle(1, 1).Face()
    right = occ.WorkPlane().MoveTo(1, 0).Rectangle(1, 1).Face()
    left.faces.name, right.faces.name = "rot", "stat"
    left.edges.name, right.edges.name = "wall_a", "wall_b"
    left.edges.Max(occ.X).name = "interface_a"
    right.edges.Min(occ.X).name = "interface_b"
    left.faces.maxh = 0.13
    right.faces.maxh = 0.21
    return Mesh(occ.OCCGeometry(occ.Compound([left, right]), dim=2).GenerateMesh(maxh=0.21))


def rotation(angle):
    c, s = np.cos(angle), np.sin(angle)
    return CF((1, 0, 0, 0, c, -s, 0, s, c), dims=(3, 3))


def rotated(mesh, angle):
    disp = GridFunction(VectorH1(mesh, order=1))
    X = CF((x, y, z))
    disp.Set(mesh.MaterialCF({"rot": rotation(angle) * X - X}, default=CF((0, 0, 0))))
    return disp


def contact_matrix(mesh, fes, integrand, disp, intorder=6, **kwargs):
    cb = ContactBoundary(mesh.Boundaries("interface_b"), mesh.Boundaries("interface_a"))
    cb.AddIntegrator(integrand)
    bf = BilinearForm(fes, check_unused=False)
    cb.Update(gf=disp, bf=bf, intorder=intorder, maxdist=0.3, quadrature="clip", **kwargs)
    bf.Assemble()
    return bf.mat


def facet_space(mesh, order):
    return Compress(FacetFESpace(mesh, order=order, definedon=mesh.Boundaries("interface_a|interface_b")))


def interior_dofs(mesh, fes, inside):
    # side-a dofs of facets away from the rim, where the two discrete rims differ
    sel = np.zeros(fes.ndof, dtype=bool)
    for el in mesh.Elements(BND):
        if mesh.GetBoundaries()[el.index] == "interface_a":
            if all(inside(mesh.ngmesh.Points()[v.nr + 1].p) for v in el.vertices):
                sel[list(fes.GetDofNrs(el))] = True
    assert sel.sum() > 10
    return sel


def check_cross_moments(mesh, fes, disp, sel, lab, moved):
    # mass of the other side and transfer of lab-frame data f, seen by side a
    # at its moved position: exact on interior facets for polynomial f
    u, v = fes.TnT()
    one = GridFunction(fes)
    one.Set(1, definedon=mesh.Boundaries("interface_a|interface_b"))
    f = GridFunction(fes)
    f.Set(lab(CF((x, y, z))), definedon=mesh.Boundaries("interface_b"))
    exact = BilinearForm(fes, check_unused=False)
    exact += u * v * ds("interface_a")
    exact.Assemble()
    m_exact = (exact.mat * one.vec).Evaluate().FV().NumPy()[sel]
    lf = LinearForm(fes)
    lf += lab(moved) * v * ds("interface_a")
    lf.Assemble()
    c_exact = lf.vec.FV().NumPy()[sel]

    A = contact_matrix(mesh, fes, u.Trace().Other() * v.Trace().Other(), disp)
    C = contact_matrix(mesh, fes, u.Trace() * v.Trace().Other(), disp)
    m = (A * one.vec).Evaluate().FV().NumPy()[sel]
    c = (C * f.vec).Evaluate().FV().NumPy()[sel]
    assert np.abs(m - m_exact).max() < 1e-12 * np.abs(m_exact).max()
    assert np.abs(c - c_exact).max() < 1e-12 * np.abs(c_exact).max()


@pytest.mark.parametrize("angle", [0.2, 0.6])
def test_clip_moments(split_cylinder, angle):
    mesh = split_cylinder
    fes = facet_space(mesh, 2)
    sel = interior_dofs(mesh, fes, lambda p: np.hypot(p[1], p[2]) < 0.65)
    check_cross_moments(mesh, fes, rotated(mesh, angle), sel,
                        lambda X: 1 + 2 * X[1] - X[2] + X[1] * X[2] - 3 * X[2] * X[2],
                        rotation(angle) * CF((x, y, z)))


@pytest.mark.parametrize("shift", [0.07, 0.31])
def test_clip_moments_2d(split_rectangle, shift):
    mesh = split_rectangle
    disp = GridFunction(VectorH1(mesh, order=1))
    disp.Set(mesh.MaterialCF({"rot": CF((0, shift))}, default=CF((0, 0))))
    fes = facet_space(mesh, 3)
    sel = interior_dofs(mesh, fes, lambda p: p[1] < 1 - shift - 0.05)
    check_cross_moments(mesh, fes, disp, sel,
                        lambda X: 1 - X[1] + 4 * X[1] ** 3,
                        CF((x, y + shift, 0)))


@pytest.mark.parametrize("angle", [0.0, 0.2, 0.6])
def test_clip_consistency(split_cylinder, angle):
    # constants are transferred exactly, and the primary weights of the pieces
    # plus the rim completion cover every primary facet
    mesh = split_cylinder
    disp = rotated(mesh, angle)
    fes = facet_space(mesh, 2)
    u, v = fes.TnT()
    one = GridFunction(fes)
    one.Set(1, definedon=mesh.Boundaries("interface_a|interface_b"))
    J = contact_matrix(mesh, fes, (u.Trace() - u.Trace().Other()) * v.Trace(), disp)
    assert Norm((J * one.vec).Evaluate()) < 1e-12
    M = contact_matrix(mesh, fes, u.Trace() * v.Trace(), disp)
    area = Integrate(1, mesh, definedon=mesh.Boundaries("interface_b"))
    assert abs(InnerProduct(M * one.vec, one.vec) - area) < 1e-12 * area


def test_clip_shared_geometry(split_cylinder):
    # a form updated with other= reuses the kept clip pieces at its own order
    mesh = split_cylinder
    disp = rotated(mesh, 0.4)
    fes = facet_space(mesh, 2)
    u, v = fes.TnT()
    mats = []
    for shared in (True, False):
        primary = ContactBoundary(mesh.Boundaries("interface_b"), mesh.Boundaries("interface_a"))
        primary.AddIntegrator(u.Trace() * v.Trace())
        primary.Update(gf=disp, bf=BilinearForm(fes, check_unused=False), intorder=2, maxdist=0.3,
                       keep_pairs=True, quadrature="clip")
        extra = {"other": primary} if shared else {}
        A = contact_matrix(mesh, fes, u.Trace().Other() * v.Trace().Other() + u.Trace() * v.Trace().Other(),
                           disp, **extra)
        mats.append(A.AsVector().FV().NumPy().copy())
    assert np.abs(mats[0] - mats[1]).max() < 1e-14 * np.abs(mats[1]).max()
