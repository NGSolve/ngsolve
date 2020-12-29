import pytest
from ngsolve import *
from netgen.geom2d import *
from netgen.csg import unit_cube

ngsglobals.msg_level=0

def test_mass_l2():
    mesh2 = Mesh(unit_square.GenerateMesh(maxh=2))
    mesh3 = Mesh(unit_cube.GenerateMesh(maxh=2))
    meshes = [mesh2, mesh3]

    for mesh in meshes:
        l2 = L2(mesh, order=3)
        c1 = FESpace([l2,l2])
        fes = FESpace([l2,FESpace([l2,l2])])

        mass = BilinearForm(fes)
        u = fes.TrialFunction()
        v = fes.TestFunction()
        mass += SymbolicBFI(u[0]*v[0] + u[1][0]*v[1][0] + u[1][1]*v[1][1])

        mass.Assemble()
        m = mass.mat

        # check if m is diagonal
        for d in range(fes.ndof):
            m[d,d] = 0.0

        assert Norm(m.AsVector())/fes.ndof<1e-15

# check if component gf keeps gf alive
def test_component_keeps_alive():
    mesh = Mesh(unit_square.GenerateMesh(maxh=2))
    assert mesh.ne == 2
    fes1 = H1(mesh)
    fes = FESpace([fes1,fes1])
    gf = GridFunction(fes)
    gf.vec[:] = 1
    gf1, gf2 = gf.components
    gf11 = gf.components[0]
    assert gf1 == gf11 # should return the same object
    assert len(gf.vec) == 8
    assert len(gf1.vec) == 4
    mesh.Refine()
    gf.Update()
    # check if gf update calls gf1 update as well
    assert len(gf.vec) == 18
    assert len(gf1.vec) == 9
    del gf
    # check if vec is still alive
    assert len(gf1.vec) == 9
    gf1.Set(0)
    for val in gf1.vec:
        assert val == 0.0

# test partially definedon number spaces
def test_number_definedon():
    geo = CSG2d()

    inner = Rectangle((0.3, 0.3), (0.4, 0.4)).Mat("inner")
    outer = Rectangle((0, 0), (1, 1))
    geo.Add((outer-inner).Mat("outer"))
    geo.Add(inner)

    mesh = Mesh(geo.GenerateMesh())

    Draw(mesh)

    h1 = H1(mesh)
    number = NumberSpace(mesh) * NumberSpace(mesh)
    number2 = NumberSpace(mesh, definedon="inner") * NumberSpace(mesh, definedon="outer")
    u = h1.TestFunction()
    t11, t12 = number.TrialFunction()
    t21, t22 = number2.TrialFunction()

    a1 = BilinearForm(testspace=h1, trialspace=number)
    a1 += u * t11 * dx("inner")
    a1 += u * t12 * dx("outer")

    a2 = BilinearForm(testspace=h1, trialspace=number2)
    a2 += u * t21 * dx("inner")
    a2 += u * t22 * dx("outer")

    a1.Assemble()
    a2.Assemble()

    for i in range(a1.mat.height):
        for j in range(a1.mat.width):
            assert a1.mat[i,j] == pytest.approx(a2.mat[i,j])

if __name__ == "__main__":
    test_component_keeps_alive()
