import numpy as np
import pathlib
import pytest
import re
from ngsolve import *
from netgen.csg import *
ngsglobals.msg_level = 0
from ngsolve.fem import LoggingCF
from ngsolve.fem import NewtonCF


def mk_fes():
    mesh = Mesh(unit_cube.GenerateMesh(maxh=1))
    return VectorH1(mesh)


def mk_cf(fes):
    return CoefficientFunction((x, y, z))


def mk_newton_cf(fes):
    _cf = mk_cf(fes)
    return NewtonCF(fes.TrialFunction() - _cf, _cf, maxiter=1)


@pytest.fixture
def fes():
    return mk_fes()


def eval_count(filename):
    with open(filename, "r") as lf:
        content = lf.read()
    evals = len(re.findall("======== Evaluate", content))
    print(evals)
    return evals


def mk_logfile(filename):
    logfile = pathlib.Path(filename)
    if logfile.exists():
        logfile.unlink(missing_ok=True)
    return logfile

# TODO: add more cases to cover all possible occurrences of caching (a lot!)


def mk_linear_form(fes, cf):
    L = LinearForm(fes)
    test_function = fes.TestFunction()
    L += InnerProduct(cf, test_function) * dx
    L += InnerProduct(cf, test_function) * ds(element_boundary=True)
    L.Assemble()
    return L


def test_cache_in_linear_form_integrator(fes):
    logfile = mk_logfile("test_cache_in_linear_form_integrator_nocache.log")
    L = mk_linear_form(fes, LoggingCF(mk_cf(fes), logfile=str(logfile)))
    baseline = eval_count(logfile)

    # covers simd version of LFIntegrator :: T_CalcElementVector (also element_vb != VOL)
    logfile = mk_logfile("test_cache_in_linear_form_integrator_cf.log")
    cL = mk_linear_form(fes, CacheCF(LoggingCF(mk_cf(fes), logfile=str(logfile))))
    assert np.allclose(cL.vec.FV().NumPy(), L.vec.FV().NumPy(), atol=1e-15, rtol=0)
    assert eval_count(logfile) == 48

    # covers non-simd version of LFIntegrator :: T_CalcElementVector (also element_vb != VOL)
    logfile = mk_logfile("test_cache_in_linear_form_integrator_ncf.log")
    cL = mk_linear_form(fes, CacheCF(LoggingCF(mk_newton_cf(fes), logfile=str(logfile))))
    assert np.allclose(cL.vec.FV().NumPy(), L.vec.FV().NumPy(), atol=1e-15, rtol=0)
    assert eval_count(logfile) == 49


def mk_bilinear_form(fes, cf):
    a = BilinearForm(fes)
    trial_function, test_function = fes.TnT()
    a += InnerProduct(cf, trial_function) * \
         InnerProduct(cf, test_function) * dx
    a += InnerProduct(cf, trial_function) * \
         InnerProduct(cf, test_function) * ds(element_boundary=True)
    a.Assemble()
    return a


def test_cache_in_bilinear_form_integrator(fes):
    logfile = mk_logfile("test_cache_in_bilinear_form_integrator_nocache.log")
    a = mk_bilinear_form(fes, LoggingCF(mk_cf(fes), logfile=str(logfile)))
    baseline = eval_count(logfile)

    # covers simd version of BLFIntegrator :: T_CalcElementMatrixAdd and T_CalcElementMatrixEBAdd
    logfile = mk_logfile("test_cache_in_bilinear_form_integrator_cf.log")
    ca = mk_bilinear_form(fes, CacheCF(LoggingCF(mk_cf(fes), logfile=str(logfile))))
    assert np.allclose(ca.mat.AsVector().FV().NumPy(), a.mat.AsVector().FV().NumPy(), atol=1e-15, rtol=0)
    assert eval_count(logfile) == 48

    # covers non-simd version of BLFIntegrator :: T_CalcElementMatrixAdd and T_CalcElementMatrixEBAdd
    logfile = mk_logfile("test_cache_in_bilinear_form_integrator_ncf.log")
    ca = mk_bilinear_form(fes, CacheCF(LoggingCF(mk_newton_cf(fes), logfile=str(logfile))))
    assert np.allclose(ca.mat.AsVector().FV().NumPy(), a.mat.AsVector().FV().NumPy(), atol=1e-15, rtol=0)
    assert eval_count(logfile) == 50


if __name__ == "__main__":
    _fes = mk_fes()
    test_cache_in_linear_form_integrator(_fes)
    test_cache_in_bilinear_form_integrator(_fes)
