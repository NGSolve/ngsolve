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
