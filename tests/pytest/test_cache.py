import numpy as np
import pathlib
import pytest
import re
from ngsolve import *
from netgen.csg import *
ngsglobals.msg_level = 0
from ngsolve.fem import LoggingCF
from ngsolve.fem import NewtonCF


# TODO: add more cases to cover all possible occurrences of caching
#  still missing:
#   * SymbolicBilinearFormIntegrator :: T_CalcElementMatrixAddShapeWise and all "Trans" versions
#   * SymbolicFacetBilinearFormIntegrator :: ApplyFromTraceValues


def mk_fes():
    mesh = Mesh(unit_cube.GenerateMesh(maxh=1))
    return VectorH1(mesh)


def mk_dg_fes():
    mesh = Mesh(unit_cube.GenerateMesh(maxh=1))
    return L2(mesh, dgjumps=True)


def mk_cf(fes):
    if fes.type == "VectorH1":
        return CoefficientFunction((x, y, z))
    else:
        return CoefficientFunction(x * y * z)


def mk_newton_cf(fes):
    _cf = mk_cf(fes)
    return NewtonCF(fes.TrialFunction() - _cf, _cf, maxiter=1)


@pytest.fixture
def fes():
    return mk_fes()


@pytest.fixture
def dg_fes():
    return mk_dg_fes()


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


def mk_linear_form(fes, cf):
    test_function = fes.TestFunction()
    L = LinearForm(fes)
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


def mk_facet_linear_form(fes, cf):
    test_function = fes.TestFunction()
    L = LinearForm(fes)
    L += InnerProduct(cf, test_function) * ds(skeleton=True)
    L.Assemble()
    return L


def test_cache_in_facet_linear_form_integrator(fes):
    logfile = mk_logfile("test_cache_in_facet_linear_form_integrator_nocache.log")
    L = mk_facet_linear_form(fes, LoggingCF(mk_cf(fes), logfile=str(logfile)))
    baseline = eval_count(logfile)

    # covers possibly existing simd version of FacetLFIntegrator :: T_CalcFacetVector
    logfile = mk_logfile("test_cache_in_facet_linear_form_integrator_cf.log")
    cL = mk_facet_linear_form(fes, CacheCF(LoggingCF(mk_cf(fes), logfile=str(logfile))))
    assert np.allclose(cL.vec.FV().NumPy(), L.vec.FV().NumPy(), atol=1e-15, rtol=0)
    assert eval_count(logfile) == 12

    # covers non-simd version of FacetLFIntegrator :: T_CalcFacetVector
    logfile = mk_logfile("test_cache_in_facet_linear_form_integrator_ncf.log")
    cL = mk_facet_linear_form(fes, CacheCF(LoggingCF(mk_newton_cf(fes), logfile=str(logfile))))
    assert np.allclose(cL.vec.FV().NumPy(), L.vec.FV().NumPy(), atol=1e-15, rtol=0)
    # NOTE: currently no simd version exists and the count for simd and non.simd are thus the same
    assert eval_count(logfile) == 12


def mk_bilinear_form(fes, cf):
    trial_function, test_function = fes.TnT()
    a = BilinearForm(fes)
    a += InnerProduct(cf, trial_function) * \
         InnerProduct(cf, test_function) * dx
    a += InnerProduct(cf, trial_function) * \
         InnerProduct(cf, test_function) * ds(element_boundary=True)
    a.Assemble()
    return a


def test_cache_in_bilinear_form_integrator(fes):
    gf = GridFunction(fes)
    gfvec = gf.vec
    resvec = gf.vec.Copy()
    gf.Interpolate(CoefficientFunction((x, y, z)))

    logfile = mk_logfile("test_cache_in_bilinear_form_integrator_nocache.log")
    a = mk_bilinear_form(fes, LoggingCF(mk_cf(fes), logfile=str(logfile)))
    a.Apply(gfvec, resvec)
    baseline = eval_count(logfile)

    # covers simd version of BLFIntegrator :: T_CalcElementMatrixAdd, T_CalcElementMatrixEBAdd and ApplyMatrix
    logfile = mk_logfile("test_cache_in_bilinear_form_integrator_cf.log")
    ca = mk_bilinear_form(fes, CacheCF(LoggingCF(mk_cf(fes), logfile=str(logfile))))
    assert eval_count(logfile) == 48
    assert np.allclose(ca.mat.AsVector().FV().NumPy(), a.mat.AsVector().FV().NumPy(), atol=1e-15, rtol=0)
    cresvec = gfvec.Copy()
    ca.Apply(gfvec, cresvec)
    assert np.allclose(cresvec.FV().NumPy(), resvec.FV().NumPy(), atol=1e-15, rtol=0)
    assert eval_count(logfile) == 96

    # covers non-simd version of BLFIntegrator :: T_CalcElementMatrixAdd, T_CalcElementMatrixEBAdd and ApplyMatrix
    logfile = mk_logfile("test_cache_in_bilinear_form_integrator_ncf.log")
    ca = mk_bilinear_form(fes, CacheCF(LoggingCF(mk_newton_cf(fes), logfile=str(logfile))))
    assert np.allclose(ca.mat.AsVector().FV().NumPy(), a.mat.AsVector().FV().NumPy(), atol=1e-15, rtol=0)
    assert eval_count(logfile) == 50
    cresvec = gfvec.Copy()
    ca.Apply(gfvec, cresvec)
    assert np.allclose(cresvec.FV().NumPy(), resvec.FV().NumPy(), atol=1e-15, rtol=0)
    assert eval_count(logfile) == 98


def mk_facet_bilinear_form(fes, cf):
    trial_function, test_function = fes.TnT()
    a = BilinearForm(fes)
    a += InnerProduct(cf, trial_function) * \
         InnerProduct(cf, test_function) * ds(skeleton=True)
    a.Assemble()
    return a


def test_cache_in_facet_bilinear_form_integrator(fes):
    gf = GridFunction(fes)
    gfvec = gf.vec
    resvec = gf.vec.Copy()
    gf.Interpolate(CoefficientFunction((x, y, z)))

    logfile = mk_logfile("test_cache_in_facet_bilinear_form_integrator_nocache.log")
    a = mk_facet_bilinear_form(fes, LoggingCF(mk_cf(fes), logfile=str(logfile)))
    a.Apply(gfvec, resvec)
    baseline = eval_count(logfile)

    # covers simd version of FacetBLFIntegrator :: T_CalcFacetMatrix, ApplyFacetMatrix (both "1-element signature")
    logfile = mk_logfile("test_cache_in_facet_bilinear_form_integrator_cf.log")
    ca = mk_facet_bilinear_form(fes, CacheCF(LoggingCF(mk_cf(fes), logfile=str(logfile))))
    assert eval_count(logfile) == 12
    assert np.allclose(ca.mat.AsVector().FV().NumPy(), a.mat.AsVector().FV().NumPy(), atol=1e-15, rtol=0)
    cresvec = gfvec.Copy()
    ca.Apply(gfvec, cresvec)
    assert np.allclose(cresvec.FV().NumPy(), resvec.FV().NumPy(), atol=1e-15, rtol=0)
    assert eval_count(logfile) == 24

    # covers non-simd version of FacetBLFIntegrator :: T_CalcFacetMatrix, ApplyFacetMatrix (both "1-element signature")
    logfile = mk_logfile("test_cache_in_facet_bilinear_form_integrator_ncf.log")
    ca = mk_facet_bilinear_form(fes, CacheCF(LoggingCF(mk_newton_cf(fes), logfile=str(logfile))))
    assert np.allclose(ca.mat.AsVector().FV().NumPy(), a.mat.AsVector().FV().NumPy(), atol=1e-15, rtol=0)
    # NOTE: currently no simd version exists and the count for simd and non.simd are thus the same
    assert eval_count(logfile) == 12
    cresvec = gfvec.Copy()
    ca.Apply(gfvec, cresvec)
    assert np.allclose(cresvec.FV().NumPy(), resvec.FV().NumPy(), atol=1e-15, rtol=0)
    assert eval_count(logfile) == 25


def mk_facet_dg_bilinear_form(dg_fes, cf):
    fes = dg_fes
    du, dv = fes.TnT()
    jump_du = du - du.Other()
    jump_dv = dv - dv.Other()

    a = BilinearForm(fes)
    a += InnerProduct(cf, cf) * jump_du * jump_dv * dx(skeleton=True)
    a.Assemble()
    return a


def test_cache_in_facet_bilinear_form_integrator_dg(dg_fes):
    fes = dg_fes
    gf = GridFunction(fes)
    gfvec = gf.vec
    resvec = gf.vec.Copy()
    gf.Interpolate(CoefficientFunction(x * y * z))

    logfile = mk_logfile("test_cache_in_facet_bilinear_form_integrator_dg_nocache.log")
    a = mk_facet_dg_bilinear_form(fes, LoggingCF(mk_cf(fes), logfile=str(logfile)))
    a.Apply(gfvec, resvec)
    baseline = eval_count(logfile)

    # covers simd version of FacetBLFIntegrator :: T_CalcFacetMatrix, ApplyFacetMatrix (both "2-element signature")
    logfile = mk_logfile("test_cache_in_facet_bilinear_form_integrator_dg_cf.log")
    ca = mk_facet_dg_bilinear_form(fes, CacheCF(LoggingCF(mk_cf(fes), logfile=str(logfile))))
    assert eval_count(logfile) == 18
    assert np.allclose(ca.mat.AsVector().FV().NumPy(), a.mat.AsVector().FV().NumPy(), atol=1e-15, rtol=0)
    cresvec = gfvec.Copy()
    ca.Apply(gfvec, cresvec)
    assert np.allclose(cresvec.FV().NumPy(), resvec.FV().NumPy(), atol=1e-15, rtol=0)
    assert eval_count(logfile) == 36

    # covers non-simd version of FacetBLFIntegrator :: T_CalcFacetMatrix, ApplyFacetMatrix (both "2-element signature")
    logfile = mk_logfile("test_cache_in_facet_bilinear_form_integrator_ncf.log")
    ca = mk_facet_dg_bilinear_form(fes, CacheCF(LoggingCF(mk_newton_cf(fes), logfile=str(logfile))))
    assert np.allclose(ca.mat.AsVector().FV().NumPy(), a.mat.AsVector().FV().NumPy(), atol=1e-15, rtol=0)
    # NOTE: currently no simd version exists and the count for simd and non.simd are thus the same
    assert eval_count(logfile) == 18
    cresvec = gfvec.Copy()
    ca.Apply(gfvec, cresvec)
    assert np.allclose(cresvec.FV().NumPy(), resvec.FV().NumPy(), atol=1e-15, rtol=0)
    assert eval_count(logfile) == 37


if __name__ == "__main__":
    _fes = mk_fes()
    test_cache_in_linear_form_integrator(_fes)
    test_cache_in_bilinear_form_integrator(_fes)
    test_cache_in_facet_linear_form_integrator(_fes)
    test_cache_in_facet_bilinear_form_integrator(_fes)

    _dg_fes = mk_dg_fes()
    test_cache_in_facet_bilinear_form_integrator_dg(_dg_fes)
