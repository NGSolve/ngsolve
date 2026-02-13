import numpy as np
import pytest
from ngsolve import *
from meshes import *


VARIANTS = [
    (False, False),
    # compile does not yet work for these edge cases
    # (True, False),
    # (True, True),
]


def _compile_fn(test_compile, real_compile):
    if not test_compile:
        return lambda cf: cf
    if real_compile:
        return lambda cf: cf.Compile(True, True)
    return lambda cf: cf.Compile()

@pytest.mark.parametrize("test_compile,real_compile", VARIANTS)
def test_gf_userdata_cache(unit_mesh_2d, test_compile, real_compile):
    """Test some edge cases for the GFCF cache in userdata"""
    mesh = unit_mesh_2d
    Compile = _compile_fn(test_compile, real_compile)

    fesr = H1(mesh, order=2)
    fesc = H1(mesh, order=2, complex=True)
    u, v = fesr.TnT()
    uc, vc = fesc.TnT()

    gfc = GridFunction(fesc)
    gfr = GridFunction(fesr)
    gfc.Set(1j)
    gfr.Set(1)

    f = LinearForm(fesr)
    f += Compile(Norm(gfc) * v) * dx

    a = BilinearForm(fesr)
    a += Compile(Norm(gfc) * u * v) * dx

    f.Assemble()
    r = f.vec.CreateVector()
    a.Apply(gfr.vec, r)
    assert np.linalg.norm(r.FV().NumPy() - f.vec.FV().NumPy()) < 1e-10

    fc = LinearForm(fesc)
    fc += Compile(gfr * vc) * dx
    fc.Assemble()

    ac = BilinearForm(fesc)
    ac += Compile(gfr * uc * vc) * dx

    rc = gfc.vec.CreateVector()
    ac.Apply(gfc.vec, rc)
    assert np.linalg.norm(rc.FV().NumPy().imag - fc.vec.FV().NumPy().real) < 1e-10

    a3 = BilinearForm(fesr)
    a3 += Compile(Norm(1j * gfr) * u * v) * dx

    r3 = gfr.vec.CreateVector()
    a3.Apply(gfr.vec, r3)
    assert np.linalg.norm(r3.FV().NumPy() - f.vec.FV().NumPy()) < 1e-10