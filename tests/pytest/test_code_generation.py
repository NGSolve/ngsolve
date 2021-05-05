from pytest import approx
from meshes import *
from ngsolve import *
ngsglobals.msg_level = 7

@pytest.mark.slow
def test_code_generation_volume_terms(unit_mesh_3d):
    fes = L2(unit_mesh_3d, order=5)
    gfu = GridFunction(fes)

    # piecewise polynomials - also test interpolation on L2 space and the resulting GridFunction
    functions = [x,y,x*y, specialcf.mesh_size, CoefficientFunction((x,y)).Norm()**2, Id(3)[:,2][2] ][-1:]

    for cf in functions:
        gfu.Set(cf)
        # should all give the same results
        cfs = [ cf.Compile(), cf.Compile(True, wait=True), gfu, gfu.Compile(), gfu.Compile(True, wait=True) ]

        for f in cfs:
            assert Integrate( (cf-f)*(cf-f), unit_mesh_3d) == approx(0)

    cf8x8 = CoefficientFunction( tuple(range(8*8)), dims=(8,8) )
    functions = [sin(x)*y, exp(x)+y*y*y, (1+x)**(1+y), cf8x8[1::3, 2:6], cf8x8[:,3], cf8x8[:,1:8:2]]
    for cf in functions:
        cfs = [ cf.Compile(), cf.Compile(True, wait=True)]
        for f in cfs:
            assert Integrate( Norm(cf-f), unit_mesh_3d) == approx(0)

    cf = atan2(1+x,1+y)
    cfs = [ cf.Compile(), cf.Compile(True, wait=True, maxderiv=0)]
    for f in cfs:
        assert Integrate( (cf-f)*(cf-f), unit_mesh_3d) == approx(0)

@pytest.mark.slow
def test_code_generation_boundary_terms(unit_mesh_3d):
    functions = [x,y,x*y, sin(x)*y, exp(x)+y*y*y, specialcf.mesh_size]
    functions = [0.1*f for f in functions]

    for cf in functions:
        # should all give the same results
        cfs = [ cf.Compile(), cf.Compile(True, wait=True)]

        for f in cfs:
            assert Integrate( (cf-f)*(cf-f), unit_mesh_3d, BND) == approx(0)

@pytest.mark.slow
def test_code_generation_volume_terms_complex(unit_mesh_3d):
    fes = L2(unit_mesh_3d, order=5, complex=True)
    gfu = GridFunction(fes)

    functions = [1J*x,1j*y,1J*x*y+y, sin(1J*x)*y+x, exp(x+1J*y)+y*y*y*sin(1J*x), 1J*specialcf.mesh_size]
    functions = functions[:1]

    for cf in functions:
        gfu.Set(cf)
        # should all give the same results
        cfs = [ cf.Compile(), cf.Compile(True, wait=True), gfu, gfu.Compile(), gfu.Compile(True, wait=True) ]

        for f in cfs:
            assert Integrate((cf-f)*Conj(cf-f), unit_mesh_3d) == approx(0)

@pytest.mark.slow
def test_code_generation_derivatives(unit_mesh_3d):
    fes = H1(unit_mesh_3d, order=4, dim=2)
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
        assert Norm(vals) == approx(0)

if __name__ == "__main__":
    test_code_generation_derivatives()
    test_code_generation_volume_terms()
    test_code_generation_volume_terms_complex()
    test_code_generation_boundary_terms()
