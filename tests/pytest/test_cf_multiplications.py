import numpy as np
import pytest
from pytest import approx
from ngsolve import *
from meshes import *


def compare_cfs(cf1, cf2, mesh):
    ad1 = np.array(cf1.dims)
    ad2 = np.array(cf2.dims)
    if ad1.size != ad2.size:
        print("'tensor-order' mismatch")
        return False
    if not np.all(ad1 == ad2):
        print("dimension mismatch")
        return False
    return Integrate(InnerProduct(cf1 - cf2, cf1 - cf2), mesh) == approx(0, abs=1e-12)


def zero_like(cf):
    if cf.dim == 1:
        return CF(0)
    else:
        return CF(tuple([0] * cf.dim), dims=tuple(cf.dims))


def test_optimizations(unit_mesh_2d):
    scal = CF(10)
    vec = CF((1, 2))
    cvec = vec.Reshape((1, 2))
    mat = CF((1, 2, 3, 4), dims=(2, 2))
    tens = CF((1, 2, 3, 4, 5, 6, 7, 8), dims=(2, 2, 2))
    I = Id(2)
    zero = CF(0)

    _compare = lambda cf1, cf2: compare_cfs(cf1, cf2, unit_mesh_2d)

    assert _compare(zero * scal, zero)
    assert _compare((zero * vec) * vec, zero)
    assert _compare(vec * (zero * vec), zero)
    assert _compare(zero * I, zero_like(I))
    assert _compare(I * zero, zero_like(I))
    assert _compare((I * zero) * I, zero_like(I))
    assert _compare(mat * (zero * vec), zero_like(vec))
    assert _compare((zero * vec).Reshape((1, 2)) * mat, zero_like(cvec))
    assert _compare((zero * cvec) * mat, zero_like(cvec))
    assert _compare(cvec * (zero * mat), zero_like(cvec))
    assert _compare((zero * cvec) * I, zero_like(cvec))
    assert _compare((mat * zero) * vec, zero_like(vec))
    assert _compare(tens * (zero * vec), zero_like(tens[:,:,0]))
    assert _compare((tens * zero) * vec, zero_like(tens[:,:,0]))

    assert _compare(I * vec, vec)
    assert _compare(cvec * I, cvec)
    assert _compare(I * mat, mat)
    assert _compare(mat * I, mat)
    assert _compare(scal * I, CF((scal, 0, 0, scal), dims=(2, 2)))
    assert _compare(I * scal, scal * I)


#@pytest.mark.skip
@pytest.mark.slow
def test_code_gen(unit_mesh_3d):
    _check_compiled = lambda cf1: compare_cfs(cf1.Compile(realcompile=True, wait=True, maxderiv=0),
                                              cf1,
                                              unit_mesh_3d)

    scal = CF(10)
    vec = CF((1, 2, 3))
    rvec = vec.Reshape((3, 1))
    mat = CF((1, 2, 3, 4, 5, 6, 7, 8, 9), dims=(3, 3))
    tens = CF(tuple(range(27)), dims=(3, 3, 3))

    assert _check_compiled(scal * vec)
    assert _check_compiled(vec * vec)
    assert _check_compiled(rvec.trans * mat)
    assert _check_compiled(mat * vec)
    assert _check_compiled(tens * vec)
    assert _check_compiled(scal * tens)
    assert np.all(np.array((scal * tens).dims) == np.array(tens.dims))


if __name__ == "__main__":
    test_optimizations(unit_mesh_2d)
    test_code_gen(unit_mesh_3d)