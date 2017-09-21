import pytest
from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from ngsolve import *
ngsglobals.msg_level = 7

def test_code_generation_volume_terms():
#     mesh = Mesh(unit_square.GenerateMesh(maxh=0.2))
    mesh = Mesh(unit_cube.GenerateMesh(maxh=0.2))
    fes = L2(mesh, order=5)
    gfu = GridFunction(fes)

    functions = [x,y,x*y, sin(x)*y, exp(x)+y*y*y, specialcf.mesh_size]

    for cf in functions:
        gfu.Set(cf)
        # should all give the same results
        cfs = [ cf.Compile(wait=True), cf.Compile(True, wait=True), gfu, gfu.Compile(wait=True), gfu.Compile(True, wait=True) ]

        for f in cfs:
            assert (Integrate ( (cf-f)*(cf-f), mesh)<1e-13)

def test_code_generation_boundary_terms():
#     mesh = Mesh(unit_square.GenerateMesh(maxh=0.2))
    mesh = Mesh(unit_cube.GenerateMesh(maxh=0.2))

    functions = [x,y,x*y, sin(x)*y, exp(x)+y*y*y, specialcf.mesh_size]
    functions = [0.1*f for f in functions]

    for cf in functions:
        # should all give the same results
        cfs = [ cf.Compile(), cf.Compile(True, wait=True)]

        for f in cfs:
            print(Integrate ( (cf-f)*(cf-f), mesh, BND))
            assert (Integrate ( (cf-f)*(cf-f), mesh, BND)<1e-13)


if __name__ == "__main__":
    test_code_generation_volume_terms()
    test_code_generation_boundary_terms()
