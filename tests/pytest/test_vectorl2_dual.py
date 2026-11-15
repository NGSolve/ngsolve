import pytest
import numpy as np
from math import log2
from ngsolve import *
from ngsolve.meshes import MakeStructured2DMesh, MakeStructured3DMesh
from netgen.geom2d import unit_square
from netgen.csg import unit_cube
import netgen.occ as occ

ngsglobals.msg_level = 0

variants = [dict(), dict(piola=True), dict(covariant=True)]
vnames = ["plain", "piola", "covariant"]


def quadmap(x, y):
    return (x + 0.2 * x * y, y + 0.1 * x * x)


def hexmap(x, y, z):
    return (x + 0.1 * y * z, y, z + 0.1 * x * y)


def mesh_trig(): return Mesh(unit_square.GenerateMesh(maxh=0.4))
def mesh_quad(): return MakeStructured2DMesh(quads=True, nx=3, ny=3)
def mesh_dquad(): return MakeStructured2DMesh(quads=True, nx=3, ny=3, mapping=quadmap)
def mesh_tet(): return Mesh(unit_cube.GenerateMesh(maxh=0.6))
def mesh_hex(): return MakeStructured3DMesh(hexes=True, nx=2, ny=2, nz=2, mapping=hexmap)
def mesh_prism(): return MakeStructured3DMesh(hexes=False, prism=True, nx=2, ny=2, nz=2)
def mesh_disk(): return Mesh(occ.OCCGeometry(occ.Circle((0, 0), 1).Face(), dim=2).GenerateMesh(maxh=0.6)).Curve(3)
def mesh_ball(): return Mesh(occ.OCCGeometry(occ.Sphere((0, 0, 0), 1)).GenerateMesh(maxh=0.8)).Curve(3)


meshes = dict(trig=mesh_trig, quad=mesh_quad, dquad=mesh_dquad, tet=mesh_tet,
              hex=mesh_hex, prism=mesh_prism, disk=mesh_disk, ball=mesh_ball)
_cache = {}


def getmesh(name):
    if name not in _cache:
        _cache[name] = meshes[name]()
    return _cache[name]


def random_gf(fes, seed=1):
    gf = GridFunction(fes)
    gf.vec.FV().NumPy()[:] = np.random.default_rng(seed).uniform(-1, 1, fes.ndof)
    return gf


@pytest.mark.parametrize("mname", list(meshes))
@pytest.mark.parametrize("vi", range(3))
@pytest.mark.parametrize("order", range(4))
def test_interpolate_reproduces(mname, vi, order):
    m = getmesh(mname)
    fes = VectorL2(m, order=order, **variants[vi])
    g0 = random_gf(fes)
    g1 = GridFunction(fes)
    g1.Interpolate(g0)
    assert np.max(np.abs(g1.vec.FV().NumPy() - g0.vec.FV().NumPy())) < 1e-12
    g2 = GridFunction(fes)
    g2.Set(g0)
    assert np.max(np.abs(g2.vec.FV().NumPy() - g1.vec.FV().NumPy())) < 1e-10
    # lazy interpolation and Set(dual=True) use the "dual" operator
    g3 = GridFunction(fes)
    g3.Set(Interpolate(g0, fes))
    assert np.max(np.abs(g3.vec.FV().NumPy() - g0.vec.FV().NumPy())) < 1e-12
    g4 = GridFunction(fes)
    g4.Set(g0, dual=True)
    assert np.max(np.abs(g4.vec.FV().NumPy() - g0.vec.FV().NumPy())) < 1e-12


@pytest.mark.parametrize("mname", ["trig", "dquad", "tet", "hex", "disk", "ball"])
@pytest.mark.parametrize("vi", [1, 2])
def test_interpolate_matches_lazy_and_dual_set(mname, vi):
    # for a general function: gf.Interpolate == Interpolate(cf, fes) == Set(dual=True)
    m = getmesh(mname)
    cf = CF((sin(x) * y, cos(y) + x)) if m.dim == 2 else CF((sin(x) * y, cos(y) + z, x * y * z))
    fes = VectorL2(m, order=2, **variants[vi])
    g1 = GridFunction(fes)
    g1.Interpolate(cf)
    g2 = GridFunction(fes)
    g2.Set(Interpolate(cf, fes))
    g3 = GridFunction(fes)
    g3.Set(cf, dual=True)
    assert np.max(np.abs(g1.vec.FV().NumPy() - g2.vec.FV().NumPy())) < 1e-12
    assert np.max(np.abs(g1.vec.FV().NumPy() - g3.vec.FV().NumPy())) < 1e-12


def dual_blocks(fes):
    u, v = fes.TnT()
    bf = BilinearForm(fes)
    bf += InnerProduct(u, v.Operator("dual")) * dx
    bf.Assemble()
    r, c, vals = bf.mat.COO()
    A = np.zeros((fes.ndof, fes.ndof))
    A[np.array(r), np.array(c)] = np.array(vals)
    return [A[np.ix_(d, d)] for d in (list(fes.GetDofNrs(el)) for el in fes.mesh.Elements(VOL))]


@pytest.mark.parametrize("mname", list(meshes))
@pytest.mark.parametrize("vi", [1, 2])
@pytest.mark.parametrize("order", range(4))
def test_dual_matrix_is_reference_mass(mname, vi, order):
    m = getmesh(mname)
    fes = VectorL2(m, order=order, **variants[vi])
    blocks = dual_blocks(fes)
    ref = blocks[0]
    assert np.max(np.abs(ref - np.diag(np.diag(ref)))) < 1e-12
    for b in blocks:
        assert np.max(np.abs(b - ref)) < 1e-12


@pytest.mark.parametrize("mname", ["trig", "dquad", "tet", "hex"])
@pytest.mark.parametrize("vi", range(3))
@pytest.mark.parametrize("order", range(3))
def test_interpolate_convergence(mname, vi, order):
    errs = []
    for n in (2, 4):
        if mname == "trig":
            m = MakeStructured2DMesh(quads=False, nx=n, ny=n, mapping=quadmap)
        elif mname == "dquad":
            m = MakeStructured2DMesh(quads=True, nx=n, ny=n, mapping=quadmap)
        elif mname == "tet":
            m = MakeStructured3DMesh(hexes=False, nx=n, ny=n, nz=n, mapping=hexmap)
        else:
            m = MakeStructured3DMesh(hexes=True, nx=n, ny=n, nz=n, mapping=hexmap)
        cf = CF((sin(2 * x) * y, cos(y) + x)) if m.dim == 2 else CF((sin(2 * x) * y, cos(y) + z, x * exp(y) * z))
        fes = VectorL2(m, order=order, **variants[vi])
        g = GridFunction(fes)
        g.Interpolate(cf)
        errs.append(sqrt(Integrate(InnerProduct(g - cf, g - cf), m, order=2 * order + 6)))
    assert log2(errs[0] / errs[1]) > order + 1 - 0.35


@pytest.mark.parametrize("vi", range(3))
def test_complex(vi):
    m = getmesh("dquad")
    fes = VectorL2(m, order=2, complex=True, **variants[vi])
    gr, gi = random_gf(VectorL2(m, order=2, **variants[vi]), 1), random_gf(VectorL2(m, order=2, **variants[vi]), 2)
    g = GridFunction(fes)
    g.Interpolate(gr + 1j * gi)
    vec = g.vec.FV().NumPy()
    assert np.max(np.abs(vec.real - gr.vec.FV().NumPy())) < 1e-12
    assert np.max(np.abs(vec.imag - gi.vec.FV().NumPy())) < 1e-12


def test_definedon():
    m = MakeStructured2DMesh(quads=False, nx=2, ny=2)
    m.ngmesh.SetMaterial(1, "a")
    fes = VectorL2(m, order=1, piola=True, definedon="a")
    g = GridFunction(fes)
    g.Interpolate(CF((1 + x, 2 - y)))
    assert sqrt(Integrate(Norm(g - CF((1 + x, 2 - y)))**2, m)) < 1e-12
