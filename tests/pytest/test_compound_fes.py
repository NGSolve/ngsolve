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

