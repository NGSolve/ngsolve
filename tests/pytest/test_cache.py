from ngsolve import *
from netgen.csg import *
ngsglobals.msg_level = 0
from ngsolve.fem import LoggingCF

import pytest


def mk_fes():
    mesh = Mesh(unit_cube.GenerateMesh(maxh=1))
    return VectorH1(mesh)


@pytest.fixture
def fes():
    return mk_fes()


def test_cache_in_linear_form_integrator(fes):
    cf = CacheCF(LoggingCF(CoefficientFunction((1,)*len(fes.components)),
                           logfile="test_cache_in_linear_form_integrator.log"))
    L = LinearForm(fes)
    L += InnerProduct(cf, CoefficientFunction(tuple(fes.TestFunction()))) * dx
    L.Assemble()


def test_cache_in_bilinear_form_integrator(fes):
    cf = CacheCF(LoggingCF(CoefficientFunction((1,)*len(fes.components)),
                           logfile="test_cache_in_bilinear_form_integrator.log"))
    a = BilinearForm(fes)
    a += InnerProduct(cf, CoefficientFunction(tuple(fes.TrialFunction()))) * \
         InnerProduct(cf, CoefficientFunction(tuple(fes.TestFunction()))) * dx
    a.Assemble()


def test_cache_in_symbolic_energy(fes):
    cf = CacheCF(LoggingCF(CoefficientFunction((1,)*len(fes.components)),
                           logfile="test_cache_in_symbolic_energy.log"))
    u = CoefficientFunction(tuple(fes.TrialFunction()))
    a = BilinearForm(fes)
    a += InnerProduct(cf, u) * InnerProduct(cf, u) * dx
    a.Assemble()


if __name__ == "__main__":
    _fes = mk_fes()
    test_cache_in_linear_form_integrator(_fes)
    test_cache_in_bilinear_form_integrator(_fes)
    test_cache_in_symbolic_energy(_fes)
