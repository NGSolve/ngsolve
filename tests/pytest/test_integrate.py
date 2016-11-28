import pytest
from ngsolve import *
from netgen.geom2d import unit_square


def test_integrate_constant():
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.2))
    cfR = []
    cfC = []
    for i in range(5):
        cfR.append(CoefficientFunction(i))
        cfC.append(CoefficientFunction(1j*i))
    cfR = CoefficientFunction(tuple(cfR))
    cfC = CoefficientFunction(tuple(cfC))
    with TaskManager():
        intR = Integrate(cfR,mesh,order=1)
        intC = Integrate(cfC,mesh,order=1)
    for i in range(5):
        assert abs(intR[i]-i) < 1e-14
        assert abs(intC[i]-1j*i) < 1e-14

def test_integrate_xy():
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.2))
    intR = Integrate(x*y,mesh)
    intC = Integrate(1j*x*y,mesh)
    assert abs(intR-1./4) < 1e-14
    assert abs(intC- 1j*1./4) < 1e-14
