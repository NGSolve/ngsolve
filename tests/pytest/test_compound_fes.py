import pytest
from ngsolve import *
from netgen.geom2d import unit_square
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

if __name__ == "__main__":
    test_component_keeps_alive()
