
import pytest
from pytest import approx
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
    
if __name__ == "__main__":
    test_pow()
    test_ParameterCF()
    test_mesh_size_cf()
    test_real()
    test_domainwise_cf()
    test_evaluate()
    test_diff()
