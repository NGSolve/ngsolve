
import pytest

@pytest.fixture
def unit_mesh_2d():
    import netgen.geom2d as g2d
    import ngsolve as ngs
    return ngs.Mesh(g2d.unit_square.GenerateMesh(maxh=0.2))

@pytest.fixture
def unit_mesh_3d():
    import netgen.csg as csg
    import ngsolve as ngs
    return ngs.Mesh(csg.unit_cube.GenerateMesh(maxh=0.2))

@pytest.fixture
def domain2_mesh_2d():
    import netgen.geom2d as g2d
    import ngsolve as ngs
    geo = g2d.SplineGeometry()
    geo.AddRectangle((-2,-2),(2,2))
    geo.AddRectangle((-1,-1),(1,1),leftdomain=2, rightdomain=1)
    geo.SetMaterial(1,"outer")
    geo.SetMaterial(2,"inner")
    return ngs.Mesh(geo.GenerateMesh(maxh=0.5))
