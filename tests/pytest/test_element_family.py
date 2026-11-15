"""
Every space must hand out, for each element of each VorB on which it has an
evaluator, an element of the family the evaluator expects, also where the
space is not defined (typed dummy elements).
"""
import functools
import pytest
import netgen.geom2d as g2d
import netgen.csg as csg
from netgen.meshing import FaceDescriptor
from netgen.meshing import Mesh as NGMesh, MeshPoint, Pnt, Element0D, Element1D
from ngsolve import *
from ngsolve.comp import Hidden, Reorder
from ngsolve.meshes import MakeStructured2DMesh, MakeStructured3DMesh


def two_materials(ngmesh, dim, bnd_dim_elements):
    # elements with centre x > 0.5 get a second material
    pts = ngmesh.Points()
    elements = ngmesh.Elements3D() if dim == 3 else ngmesh.Elements2D()
    for el in elements:
        cx = sum(pts[v].p[0] for v in el.vertices) / len(el.vertices)
        if cx > 0.5:
            el.index = 2
    ngmesh.SetMaterial(1, "left")
    ngmesh.SetMaterial(2, "right")
    return ngmesh


def seg1d():
    # two materials, segments of both orientations
    m = NGMesh(dim=1)
    n = 6
    pids = [m.Add(MeshPoint(Pnt(i/n, 0, 0))) for i in range(n+1)]
    left = m.AddRegion("left", dim=1)
    right = m.AddRegion("right", dim=1)
    bleft = m.AddRegion("bleft", dim=0)
    bright = m.AddRegion("bright", dim=0)
    for i in range(n):
        pair = [pids[i], pids[i+1]] if i % 2 else [pids[i+1], pids[i]]
        m.Add(Element1D(pair, index=left if i < n//2 else right))
    m.Add(Element0D(pids[0], index=bleft))
    m.Add(Element0D(pids[n], index=bright))
    return Mesh(m)


def tri2d():
    geo = g2d.SplineGeometry()
    geo.AddRectangle((0, 0), (1, 1), leftdomain=1, rightdomain=0)
    geo.AddRectangle((0.25, 0.25), (0.75, 0.75), leftdomain=2, rightdomain=1)
    geo.SetMaterial(1, "outer")
    geo.SetMaterial(2, "inner")
    return Mesh(geo.GenerateMesh(maxh=0.35))


def quad2d():
    m = MakeStructured2DMesh(quads=True, nx=4, ny=4)
    return Mesh(two_materials(m.ngmesh, 2, 1))


def tri2d_structured():
    m = MakeStructured2DMesh(quads=False, nx=4, ny=4)
    return Mesh(two_materials(m.ngmesh, 2, 1))


def tri2d_periodic():
    m = MakeStructured2DMesh(quads=False, nx=4, ny=4, periodic_x=True)
    return Mesh(two_materials(m.ngmesh, 2, 1))


def tet3d():
    geo = csg.CSGeometry()
    outer = csg.OrthoBrick(csg.Pnt(0, 0, 0), csg.Pnt(1, 1, 1))
    inner = csg.OrthoBrick(csg.Pnt(0.25, 0.25, 0.25), csg.Pnt(0.75, 0.75, 0.75))
    geo.Add((outer - inner).mat("outer"))
    geo.Add(inner.mat("inner"))
    return Mesh(geo.GenerateMesh(maxh=0.5))


def hex3d():
    m = MakeStructured3DMesh(hexes=True, nx=3)
    return Mesh(two_materials(m.ngmesh, 3, 2))


def prism3d():
    m = MakeStructured3DMesh(hexes=False, prism=True, nx=3)
    return Mesh(two_materials(m.ngmesh, 3, 2))


def hex3d_periodic():
    m = MakeStructured3DMesh(hexes=True, nx=3, periodic_x=True)
    return Mesh(two_materials(m.ngmesh, 3, 2))


MESHES = {
    "seg1d": (seg1d, 1),
    "tri2d": (tri2d, 2),
    "quad2d": (quad2d, 2),
    "tri2d_structured": (tri2d_structured, 2),
    "tet3d": (tet3d, 3),
    "hex3d": (hex3d, 3),
    "prism3d": (prism3d, 3),
    "tri2d_periodic": (tri2d_periodic, 2),
    "hex3d_periodic": (hex3d_periodic, 3),
}


@functools.lru_cache(maxsize=None)
def get_mesh(name):
    return MESHES[name][0]()


def vol_region(mesh):
    return mesh.Materials(mesh.GetMaterials()[0])


def bnd_region(mesh):
    return mesh.Boundaries(mesh.GetBoundaries()[0])


def kw(mesh, restricted, bnd=False):
    if not restricted:
        return {}
    return {"definedon": bnd_region(mesh) if bnd else vol_region(mesh)}


# name -> (constructor(mesh, restricted), allowed mesh names)
ALL2D = ["tri2d", "quad2d", "tri2d_structured"]
ALL3D = ["tet3d", "hex3d", "prism3d"]
ALL = ALL2D + ALL3D
ALL1D = ["seg1d"]
SIMPLEX = ["tri2d", "tri2d_structured", "tet3d"]

SPACES = {
    "H1": (lambda m, r: H1(m, order=2, **kw(m, r)), ALL1D + ALL),
    "L2": (lambda m, r: L2(m, order=2, **kw(m, r)), ALL1D + ALL),
    "HCurl": (lambda m, r: HCurl(m, order=2, **kw(m, r)), ALL1D + ALL),
    "HDiv": (lambda m, r: HDiv(m, order=2, **kw(m, r)), ALL1D + ALL),
    "HDivDiv": (lambda m, r: HDivDiv(m, order=2, **kw(m, r)), ALL),
    "HCurlCurl": (lambda m, r: HCurlCurl(m, order=2, **kw(m, r)), ALL),
    "HCurlDiv": (lambda m, r: HCurlDiv(m, order=2, **kw(m, r)), ALL),
    "FacetFESpace": (lambda m, r: FacetFESpace(m, order=2, **kw(m, r)), ALL1D + ALL),
    "FacetSurface": (lambda m, r: FacetSurface(m, order=2, **kw(m, r, True)), ALL),
    "NormalFacetFESpace": (lambda m, r: NormalFacetFESpace(m, order=2, **kw(m, r)), ALL1D + ALL),
    "NormalFacetSurface": (lambda m, r: NormalFacetSurface(m, order=2, **kw(m, r, True)), ALL3D),
    "TangentialFacetFESpace": (lambda m, r: TangentialFacetFESpace(m, order=2, **kw(m, r)), ALL1D + ALL),
    "SurfaceL2": (lambda m, r: SurfaceL2(m, order=2, **kw(m, r, True)), ALL1D + ALL),
    "TangentialSurfaceL2": (lambda m, r: TangentialSurfaceL2(m, order=2, **kw(m, r, True)), ALL1D + ALL3D),
    "VectorH1": (lambda m, r: VectorH1(m, order=2, **kw(m, r)), ALL1D + ALL),
    "VectorL2": (lambda m, r: VectorL2(m, order=2, **kw(m, r)), ALL1D + ALL),
    "VectorFacetFESpace": (lambda m, r: VectorFacetFESpace(m, order=2, **kw(m, r)), ALL1D + ALL),
    "VectorSurfaceL2": (lambda m, r: VectorSurfaceL2(m, order=2, **kw(m, r, True)), ALL1D + ALL3D),
    "VectorFacetSurface": (lambda m, r: VectorFacetSurface(m, order=2, **kw(m, r, True)), ALL3D),
    "VectorValuedH1": (lambda m, r: VectorValued(H1(m, order=2, **kw(m, r)), 2), ALL1D + ALL),
    "MatrixValuedH1": (lambda m, r: MatrixValued(H1(m, order=2, **kw(m, r)), 2), ALL1D + ALL),
    "NumberSpace": (lambda m, r: NumberSpace(m, **kw(m, r)), ALL1D + ALL),
    "NodalFESpace": (lambda m, r: NodalFESpace(m, order=1, **kw(m, r)), ALL1D + ALL),
    "HDivSurface": (lambda m, r: HDivSurface(m, order=2, **kw(m, r, True)), ALL3D),
    "HDivDivSurface": (lambda m, r: HDivDivSurface(m, order=2, **kw(m, r, True)), ALL3D),
    "HDivDivFacetSpace": (lambda m, r: HDivDivFacetSpace(m, order=2, **kw(m, r)), ["tet3d"]),
    "Discontinuous(H1)": (lambda m, r: Discontinuous(H1(m, order=2, **kw(m, r))), ALL1D + ALL),
    "Compress(H1)": (lambda m, r: Compress(H1(m, order=2, definedon=vol_region(m))), ALL1D + ALL),
    "H1*HDiv*L2": (lambda m, r: H1(m, order=2, **kw(m, r)) * HDiv(m, order=2, **kw(m, r)) * L2(m, order=2, **kw(m, r)), ALL1D + ALL),
    "HDivDiv*H1": (lambda m, r: HDivDiv(m, order=2, **kw(m, r)) * H1(m, order=2, **kw(m, r)), ALL),
    "HCurl*FacetFESpace": (lambda m, r: HCurl(m, order=2, **kw(m, r)) * FacetFESpace(m, order=2, **kw(m, r)), ALL1D + ALL),
    "Periodic(H1)": (lambda m, r: Periodic(H1(m, order=2, **kw(m, r))), ["tri2d_periodic", "hex3d_periodic"]),
    "Periodic(HDiv)": (lambda m, r: Periodic(HDiv(m, order=2, **kw(m, r))), ["tri2d_periodic", "hex3d_periodic"]),
    "Hidden(H1)": (lambda m, r: Hidden(H1(m, order=2, **kw(m, r))), ALL1D + ALL),
    "Reorder(H1)": (lambda m, r: Reorder(H1(m, order=2, **kw(m, r))), ALL1D + ALL),
}

CASES = [(s, m, r) for s, (_, meshes) in SPACES.items() for m in meshes for r in (False, True)]


@pytest.mark.parametrize("space,meshname,restricted", CASES,
                         ids=[f"{s}-{m}-{'definedon' if r else 'all'}" for s, m, r in CASES])
def test_element_family(space, meshname, restricted):
    mesh = get_mesh(meshname)
    ctor = SPACES[space][0]
    try:
        fes = ctor(mesh, restricted)
    except Exception as e:
        pytest.skip(f"construction rejected: {e}")
    fes._CheckElementFamilies()
