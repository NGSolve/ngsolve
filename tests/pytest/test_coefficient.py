import pytest
from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from ngsolve import *
from netgen.geom2d import SplineGeometry



def test_ParameterCF():
    p = Parameter(23)
    assert type(p) == Parameter
                    
def test_mesh_size_cf():
    for mesh in [ Mesh(unit_cube.GenerateMesh(maxh=0.2)), Mesh(unit_square.GenerateMesh(maxh=0.2))]:
        dim = mesh.dim
        fes = L2(mesh, order=0)
        gfu = GridFunction(fes)
        cf = specialcf.mesh_size 
        if dim==2:
            gfu.Set(cf*cf/2)
        else:
            gfu.Set(cf*cf*cf/6)
        assert abs( sum(gfu.vec)-1 ) < 1e-12

        bfi = BilinearForm(fes)
        bfi += SymbolicBFI( fes.TrialFunction()*(specialcf.mesh_size - specialcf.mesh_size.Compile(True, wait=True))*fes.TestFunction(), element_boundary=True)
        bfi.Assemble()
        v = bfi.mat.AsVector().FV()
        for val in v:
            assert abs(val) < 1e-14

def test_real():
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.2))
    cf = CoefficientFunction(1+2j)
    assert cf.real(mesh(0.4,0.4)) == 1
    assert cf.imag(mesh(0.2,0.6)) == 2

def test_domainwise_cf():
    geo = SplineGeometry()
    geo.AddCircle ( (0, 0), r=1, leftdomain=1, rightdomain=2)
    geo.AddCircle ( (0, 0), r=1.4, leftdomain=2, rightdomain=0)
    mesh = Mesh(geo.GenerateMesh(maxh=0.1))
    c_vec = CoefficientFunction([(x,y),(1,3)])
    c = c_vec[0]*c_vec[1]

    c_false = c.Compile(False);
    error_false = Integrate((c-c_false)*(c-c_false), mesh)
    assert abs(error_false) < 1e-14

    c_true = c.Compile(True, wait=True);
    error_true = Integrate((c-c_true)*(c-c_true), mesh)
    assert abs(error_true) < 1e-14

if __name__ == "__main__":
    test_ParameterCF()
    test_mesh_size_cf()
    test_real()
    test_domainwise_cf()
