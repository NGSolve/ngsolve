import pytest
from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from ngsolve import *
ngsglobals.msg_level = 7

def test_code_generation_volume_terms():
    mesh = Mesh(unit_cube.GenerateMesh(maxh=0.2))
    fes = L2(mesh, order=5)
    gfu = GridFunction(fes)

    functions = [x,y,x*y, sin(x)*y, exp(x)+y*y*y, specialcf.mesh_size, CoefficientFunction((x,y)).Norm()**2]

    for cf in functions:
        gfu.Set(cf)
        # should all give the same results
        cfs = [ cf.Compile(), cf.Compile(True, wait=True), gfu, gfu.Compile(), gfu.Compile(True, wait=True) ]

        for f in cfs:
            assert (Integrate ( (cf-f)*(cf-f), mesh)<1e-13)

def test_code_generation_boundary_terms():
    mesh = Mesh(unit_cube.GenerateMesh(maxh=0.2))

    functions = [x,y,x*y, sin(x)*y, exp(x)+y*y*y, specialcf.mesh_size]
    functions = [0.1*f for f in functions]

    for cf in functions:
        # should all give the same results
        cfs = [ cf.Compile(), cf.Compile(True, wait=True)]

        for f in cfs:
            print(Integrate ( (cf-f)*(cf-f), mesh, BND))
            assert (Integrate ( (cf-f)*(cf-f), mesh, BND)<1e-13)

def test_code_generation_volume_terms_complex():
    mesh = Mesh(unit_cube.GenerateMesh(maxh=0.2))
    fes = L2(mesh, order=5, complex=True)
    gfu = GridFunction(fes)

    functions = [1J*x,1j*y,1J*x*y+y, sin(1J*x)*y+x, exp(x+1J*y)+y*y*y*sin(1J*x), 1J*specialcf.mesh_size]
    functions = functions[:1]

    for cf in functions:
        gfu.Set(cf)
        print ("norm gfu.vec = ", Norm(gfu.vec))
        # should all give the same results
        cfs = [ cf.Compile(), cf.Compile(True, wait=True), gfu, gfu.Compile(), gfu.Compile(True, wait=True) ]

        for f in cfs:
            print ("integral f =", Integrate(f, mesh))
            print ("integral cf =", Integrate(cf, mesh))
            print ("integral cf-f =", Integrate(cf-f, mesh))
            print ("integral Conj(cf-f) =", Integrate(Conj(cf-f), mesh))
            error = (cf-f)*Conj(cf-f)
            assert (abs(Integrate ( error, mesh))<1e-13)

def test_code_generation_derivatives():
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.2))
    fes = H1(mesh, order=4, dim=2)
    gfu = GridFunction(fes)
    u,v = fes.TnT()

    gfu.Set(CoefficientFunction((x*x+y*y*3, x*y)))
    cf = Norm(u[0]*u) + Norm(u[0]*grad(u))
    cfs = [cf.Compile(), cf.Compile(True, wait=True)]

    aref = BilinearForm(fes, symmetric=False)
    aref += SymbolicEnergy( cf )
    aref.AssembleLinearization(gfu.vec)
    vals_ref = aref.mat.AsVector()

    for f in cfs:
        a = BilinearForm(fes)
        a += SymbolicEnergy( f )
        a.AssembleLinearization(gfu.vec)
        vals = a.mat.AsVector()
        vals -= vals_ref
        assert Norm(vals) < 1e-13

if __name__ == "__main__":
    test_code_generation_derivatives()
    test_code_generation_volume_terms()
    test_code_generation_volume_terms_complex()
    test_code_generation_boundary_terms()
