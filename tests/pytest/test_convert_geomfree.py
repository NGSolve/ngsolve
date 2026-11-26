import pytest
from ngsolve import *
from netgen.occ import *


def meshes():
    box = Box((0,0,0), (1,1,1))
    yield Mesh(OCCGeometry(box).GenerateMesh(maxh=0.5))
    yield Mesh(OCCGeometry(box - Sphere((0.5,0.5,0.5), 0.25)).GenerateMesh(maxh=0.4)).Curve(3)


pairs = [
    lambda m, k: (L2(m, order=k), FacetFESpace(m, order=k)),
    lambda m, k: (HCurl(m, order=k), TangentialFacetFESpace(m, order=k)),
    lambda m, k: (VectorL2(m, order=k, covariant=True), TangentialFacetFESpace(m, order=k)),
    lambda m, k: (HDiv(m, order=k), NormalFacetFESpace(m, order=k)),
    lambda m, k: (VectorL2(m, order=k, piola=True), NormalFacetFESpace(m, order=k)),
    lambda m, k: (HDivDiv(m, order=k), Discontinuous(HDivDiv(m, order=k))),
]


@pytest.mark.parametrize("pair", pairs)
@pytest.mark.parametrize("order", [1, 2])
def test_geom_free_equals_sparse(pair, order):
    for mesh in meshes():
        fa, fb = pair(mesh, order)
        u = GridFunction(fa).vec
        u.SetRandom(1)
        v1 = GridFunction(fb).vec.CreateVector()
        v2 = v1.CreateVector()
        v1.data = ConvertOperator(fa, fb) * u
        v2.data = ConvertOperator(fa, fb, geom_free=True) * u
        assert Norm(v1) > 0
        assert Norm(v1-v2) < 1e-12 * Norm(v1)


@pytest.mark.parametrize("order", [0, 1, 2, 3])
def test_normal_trace_tet(order):
    n = specialcf.normal(3)
    for mesh in meshes():
        for fa in [HDiv(mesh, order=order), HDiv(mesh, order=order, RT=True)]:
            u = GridFunction(fa)
            u.vec.SetRandom(1)
            uhat = GridFunction(NormalFacetFESpace(mesh, order=order))
            uhat.vec.data = ConvertOperator(fa, uhat.space) * u.vec
            err = Integrate(((u-uhat)*n)**2 * dx(element_boundary=True), mesh)
            assert err < 1e-20 * Integrate((u*n)**2 * dx(element_boundary=True), mesh)


def test_tangential_trace_tet():
    n = specialcf.normal(3)
    tang = lambda v: v - (v*n)*n
    for mesh in meshes():
        fa = HCurl(mesh, order=2)
        u = GridFunction(fa)
        u.vec.SetRandom(1)
        uhat = GridFunction(TangentialFacetFESpace(mesh, order=2))
        uhat.vec.data = ConvertOperator(fa, uhat.space, geom_free=True) * u.vec
        err = Integrate(tang(u-uhat)**2 * dx(element_boundary=True), mesh)
        assert err < 1e-20 * Integrate(tang(u)**2 * dx(element_boundary=True), mesh)
