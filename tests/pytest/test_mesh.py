from ngsolve import *
from netgen.csg import *
ngsglobals.msg_level = 0

def test_multiple_meshes_refine():
    mesh = Mesh(unit_cube.GenerateMesh(maxh=1))
    mesh.Refine()
    parents = [mesh.GetParentVertices(v) for v in range(mesh.nv)]
    print(parents)

    mesh2 = Mesh(unit_cube.GenerateMesh(maxh=0.5))
    mesh2.Refine()

    for v in range(mesh.nv):
        assert parents[v] == mesh.GetParentVertices(v)

def test_call_operator_with_vectorization():
    mesh = Mesh(unit_cube.GenerateMesh(maxh=1))
    p = mesh(0.5,0.5,0.5)
    p2 = mesh([0.5, 0.1],0.5,0.5)

def test_neighbours():
    geo = CSGeometry()

    bot = Plane(Pnt(0,0,0), Vec(0,0,-1)).bc("bot")
    top = Plane(Pnt(10,10,10), Vec(0,0,1)).bc("top")
    box = OrthoBrick(Pnt(0,0,-1), Pnt(10,10,11)).bc("rest") * top * bot

    plate = OrthoBrick(Pnt(3,3,4.8), Pnt(7,7,5.2))
    cutx = Plane(Pnt(5,5,5), Vec(1,0,0)).bc("cutx")
    cuty = Plane(Pnt(5,5,5), Vec(0,1,0)).bc("cuty")

    cond1 = plate * cutx * cuty
    cond2 = plate * cutx - cuty
    cond3 = plate * cuty - cutx
    cond4 = plate - cutx - cuty

    geo.Add(cond1.mat("cond1"))
    geo.Add(cond2.mat("cond2"))
    geo.Add(cond3.mat("cond3"))
    geo.Add(cond4.mat("cond4"))

    geo.Add((box-plate).mat("air"))

    mesh = Mesh(geo.GenerateMesh(maxh=2))

    xpl = mesh.Boundaries("cutx")
    nmats = [mesh.GetMaterials()[i] for i, val in enumerate(xpl.Neighbours(VOL).Mask()) if val]
    print(nmats)
    assert nmats == [ "cond1", "cond2", "cond3", "cond4"]
    cond12 = mesh.Materials("cond1")
    print(set(mesh.GetBoundaries()[i] for i, val in enumerate(cond12.Boundaries().Mask()) if val))
    assert set((mesh.GetBoundaries()[i] for i, val in enumerate(cond12.Boundaries().Mask()) if val)) == {"default", "cuty", "cutx"}
    print(set(mesh.GetBoundaries()[i] for i, val in enumerate(mesh.Materials("air").Boundaries().Mask()) if val))
    assert not any([bnd in set(mesh.GetBoundaries()[i] for i, val in enumerate(mesh.Materials("air").Boundaries().Mask()) if val) for bnd in ["cuty", "cutx"]])

    import numpy
    midedge = numpy.prod([reg.Neighbours(BBND) for reg in mesh.Materials("cond.*").Split()])
    print(midedge.Mask())
    assert sum(midedge.Mask()) == 1

    midedge2 = mesh.Boundaries("cuty").Boundaries() * mesh.Boundaries("cutx").Boundaries()
    print(midedge2.Mask())
    assert midedge == midedge2
    print(mesh.Materials("cond1").Neighbours(VOL).Mask())
    assert mesh.Materials("cond1").Neighbours(VOL) == mesh.Materials("cond2|cond3|air")

if __name__ == "__main__":
    test_neighbours()
