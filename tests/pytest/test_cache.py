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
        lines = lf.readlines()
    lines = [l for l in lines if l.startswith("======== Evaluate")]
    evals = 0
    
    # have_simd = len(re.findall("======== Evaluate\\(ngfem::SIMD_Mapped", lines))
    have_simd = False
    for line in lines:
        if "::Mapped" in line:
            evals = 1 if have_simd else evals+1
            have_simd = False
        else:
            have_simd = True
            evals += 1
    # evals_regular = len(re.findall("======== Evaluate\\(ngfem::Mapped", content))
    # evals_simd = len(re.findall("======== Evaluate\\(ngfem::SIMD_Mapped", content))
    # evals = max(evals_simd, evals_regular)
    print(evals)
    return evals


_logfile_counter = 0
def mk_logfile(filename):
    global _logfile_counter
    logfile = pathlib.Path(f"{filename}_{_logfile_counter}")
    if logfile.exists():
        logfile.unlink()
    _logfile_counter += 1
    return logfile


def mk_linear_form(fes, cf):
    test_function = fes.TestFunction()
    L = LinearForm(fes)
    L += InnerProduct(cf, test_function) * dx
    L += InnerProduct(cf, test_function) * ds(element_boundary=True)
    L.Assemble()
    return L

def compare_lf(make_lf, fes):
    logfile = mk_logfile("baseline.log")
    L = make_lf(fes, LoggingCF(mk_cf(fes), logfile=str(logfile)))
    baseline = eval_count(logfile)

    for cf in [mk_cf(fes), mk_newton_cf(fes)]:
        logfile = mk_logfile("check.log")
        cL = make_lf(fes, CacheCF(LoggingCF(cf, logfile=str(logfile))))
        assert np.allclose(cL.vec.FV().NumPy(), L.vec.FV().NumPy(), atol=1e-15, rtol=0)
        assert eval_count(logfile) < baseline


def test_cache_in_linear_form_integrator(fes):
    compare_lf(mk_linear_form, fes)


def mk_facet_linear_form(fes, cf):
    test_function = fes.TestFunction()
    L = LinearForm(fes)
    L += InnerProduct(cf, test_function) * ds(skeleton=True)
    L.Assemble()
    return L


def test_cache_in_facet_linear_form_integrator(fes):
    compare_lf(mk_facet_linear_form, fes)


def mk_bilinear_form(fes, cf):
    trial_function, test_function = fes.TnT()
    a = BilinearForm(fes)
    a += InnerProduct(cf, trial_function) * \
         InnerProduct(cf, test_function) * dx
    a += InnerProduct(cf, trial_function) * \
         InnerProduct(cf, test_function) * ds(element_boundary=True)
    a.Assemble()
    return a


def compare_blf(make_blf, fes):
    logfile = mk_logfile("baseline.log")
    a = make_blf(fes, LoggingCF(mk_cf(fes), logfile=str(logfile)))
    baseline = eval_count(logfile)

    gf = GridFunction(fes)
    gfvec = gf.vec
    resvec = gf.vec.Copy()
    # gf.Interpolate(CoefficientFunction((x, y, z)))

    for cf in [mk_cf(fes), mk_newton_cf(fes)]:
        logfile = mk_logfile("check.log")
        ca = make_blf(fes, CacheCF(LoggingCF(cf, logfile=str(logfile))))
        count_assemble = eval_count(logfile)
        assert count_assemble < baseline
        assert np.allclose(ca.mat.AsVector().FV().NumPy(), a.mat.AsVector().FV().NumPy(), atol=1e-15, rtol=0)
        cresvec = gfvec.Copy()
        ca.Apply(gfvec, cresvec)
        assert np.allclose(cresvec.FV().NumPy(), resvec.FV().NumPy(), atol=1e-15, rtol=0)
        assert eval_count(logfile)-count_assemble < baseline


def test_cache_in_bilinear_form_integrator(fes):
    compare_blf(mk_bilinear_form, fes)

def mk_facet_bilinear_form(fes, cf):
    trial_function, test_function = fes.TnT()
    a = BilinearForm(fes)
    a += InnerProduct(cf, trial_function) * \
         InnerProduct(cf, test_function) * ds(skeleton=True)
    a.Assemble()
    return a


def test_cache_in_facet_bilinear_form_integrator(fes):
    compare_blf(mk_facet_bilinear_form, fes)


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
    compare_blf(mk_facet_dg_bilinear_form, dg_fes)


if __name__ == "__main__":
    _fes = mk_fes()
    test_cache_in_linear_form_integrator(_fes)
    test_cache_in_bilinear_form_integrator(_fes)
    test_cache_in_facet_linear_form_integrator(_fes)
    test_cache_in_facet_bilinear_form_integrator(_fes)

    _dg_fes = mk_dg_fes()
    test_cache_in_facet_bilinear_form_integrator_dg(_dg_fes)
