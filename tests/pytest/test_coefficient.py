
import pytest
from pytest import approx
from math import isnan
from ngsolve import *
from meshes import *

def test_ParameterCF():
    p = Parameter(23)
    assert type(p) == Parameter
                    
def test_mesh_size_cf(unit_mesh_3d, unit_mesh_2d):
    for mesh in [ unit_mesh_3d, unit_mesh_2d]:
        dim = mesh.dim
        fes = L2(mesh, order=0)
        gfu = GridFunction(fes)
        cf = specialcf.mesh_size 
        if dim==2:
            gfu.Set(cf*cf/2)
        else:
            gfu.Set(cf*cf*cf/6)
        assert sum(gfu.vec) == approx(1,rel=1e-12)

        bfi = BilinearForm(fes)
        bfi += SymbolicBFI( fes.TrialFunction()*(specialcf.mesh_size - specialcf.mesh_size.Compile(True, wait=True))*fes.TestFunction(), element_boundary=True)
        bfi.Assemble()
        v = bfi.mat.AsVector().FV().NumPy()
        assert all((val == approx(0) for val in v))

def test_real(unit_mesh_2d):
    cf = CoefficientFunction(1+2j)
    assert cf.real(unit_mesh_2d(0.4,0.4)) == 1
    assert cf.imag(unit_mesh_2d(0.2,0.6)) == 2

def test_pow(unit_mesh_2d):
    base = (x+0.1)
    CompareCfs = lambda c1, c2, mesh: Integrate((c1-c2)*(c1-c2), mesh) == approx(0,abs=1e-12)
    for p in range(10):
        c = CoefficientFunction(1.0)
        for i in range(p):
            c = c*base
        assert CompareCfs(base**p,c,unit_mesh_2d)
        assert CompareCfs(base**(-p),1.0/c,unit_mesh_2d)
        assert CompareCfs(base**float(p),c,unit_mesh_2d)
        assert CompareCfs(base**CoefficientFunction(p),c,unit_mesh_2d)

def test_pow_negative_base(unit_mesh_2d):
    # float exponent on a sign-changing base
    assert Integrate((0.5-x)**2., unit_mesh_2d) == approx(1/12, rel=1e-12)
    assert Integrate(((0.5-x)**2.).Compile(), unit_mesh_2d) == approx(1/12, rel=1e-12)
    base = 0.5-x
    CompareCfs = lambda c1, c2, mesh: Integrate((c1-c2)*(c1-c2), mesh) == approx(0,abs=1e-12)
    for p in range(1,6):
        c = CoefficientFunction(1.0)
        for i in range(p):
            c = c*base
        assert CompareCfs(base**float(p),c,unit_mesh_2d)
        assert CompareCfs(base**CoefficientFunction(float(p)),c,unit_mesh_2d)
        assert CompareCfs((base**float(p)).Compile(),c,unit_mesh_2d)
        # Parameter exponent is not reduced to products -> SIMD pow
        assert CompareCfs(base**Parameter(float(p)),c,unit_mesh_2d)
    assert Integrate(base**Parameter(0.), unit_mesh_2d) == approx(1)
    assert isnan(Integrate(base**Parameter(1.5), unit_mesh_2d))

def test_domainwise_cf(domain2_mesh_2d):
    c_vec = CoefficientFunction([(x,y),(1,3)])
    c = c_vec[0]*c_vec[1]

    assert Integrate(c_vec, domain2_mesh_2d, definedon=domain2_mesh_2d.Materials("inner")) == approx([4,12])

    c_false = c.Compile(False);
    error_false = Integrate((c-c_false)*(c-c_false), domain2_mesh_2d)
    assert error_false == approx(0)

    c_true = c.Compile(True, wait=True);
    error_true = Integrate((c-c_true)*(c-c_true), domain2_mesh_2d)
    assert error_true == approx(0)

def test_evaluate(unit_mesh_2d):
    import numpy as np
    pnts = np.linspace(0.1,0.9,9)
    cf = CoefficientFunction((x,y))
    cf2 = CoefficientFunction((y, x * 1J))
    mips = unit_mesh_2d(pnts,0.5)
    vals = cf(mips)
    vals2 = cf2(mips)
    assert vals == approx(np.array(list(zip(pnts,[0.5]*10))))
    assert vals2 == approx(np.array(list(zip([0.5 + 0J] * 10, pnts*1J))))
    assert x(unit_mesh_2d(0.5,0.5)) == approx(0.5)


def test_diff(unit_mesh_3d):
    u = CF( (x ** 2 * z * y-0.08, y ** 3 * x, z ** 2 * x * y) ).MakeVariable()
    F = CF( (exp(x+y),sin(x-z),cos(y**2+z*x), x**2-y,log( 2+0.2*x*y),sqrt(4+y+z), u), dims=(3,3) ) + 10*Id(3)

    cfs = [Det(F), Cof(F), Sym(F), Skew(F), Inv(F), F.trans*F, F*u, InnerProduct(F,F), InnerProduct(u,F*u), IfPos(u[0], F, Inv(F)), Norm(F), Inv(F)-F+ Cof(F), Cross(u,F*u), 1/(Norm(u)+10)*F, F[:,0],F.Compile(), F.ExtendDimension(dims=(4,3))*u]

    for cf in cfs:
        assert Integrate( Norm(cf.Diff(u,CF((1,0,0)))-cf.Diff(u)*CF((1,0,0))),unit_mesh_3d) == approx(0.0)

def test_diff_nested(unit_mesh_2d):
    # Jacobian of Jacobian vs directional derivatives, matrix and vector variables
    G = CF(tuple(0.1*i+x*(i%3)+y*(i%2)+(i in (0,4,8)) for i in range(9)), dims=(3,3)).MakeVariable()
    v = CF((x+1,y,2)).MakeVariable()
    B = CF((1,x,2, y,3,1, x,0,4), dims=(3,3))
    def unit(var, k):
        return CF(tuple(float(i==k) for i in range(var.dim)), dims=tuple(var.dims))
    cases = [ (Cof(G), G, G), (G*Cof(G), G, G), (G*G*v, G, v), (Inv(G)*B, G, G),
              (Cross(v,G*v), v, v), (Cross(G*v,v), G, v), (Cross(Cof(G)*v,G*v), G, G),
              (OuterProduct(G*v,v), v, G) ]
    mp = unit_mesh_2d(0.3, 0.4)
    for cf, x1, x2 in cases:
        H = cf.Diff(x1).Diff(x2)
        h = H(mp)
        for k in range(x1.dim):
            for l in range(x2.dim):
                ref = cf.Diff(x1,unit(x1,k)).Diff(x2,unit(x2,l))(mp)
                ref = ref if isinstance(ref, tuple) else (ref,)
                hkl = h[(k*x2.dim+l)::x1.dim*x2.dim]
                assert hkl == approx(ref, abs=1e-12)

if __name__ == "__main__":
    test_pow()
    test_ParameterCF()
    test_mesh_size_cf()
    test_real()
    test_domainwise_cf()
    test_evaluate()
    test_diff()
