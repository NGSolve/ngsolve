import pytest
import numpy as np

import ngsolve.bem as bem
import ngsolve.fmm as fmm
from ngsolve.fmm import FMM, SingularMLExpansion

rng = np.random.default_rng(20261112)
N, M = 300, 200
SRC = rng.random((N, 3))
TGT = rng.random((M, 3)) + [0.3, 0, 0]
NRM = rng.standard_normal((N, 3))
NRM /= np.linalg.norm(NRM, axis=1)[:, None]
Q = rng.standard_normal(N) + 1j * rng.standard_normal(N)
Q3 = rng.standard_normal((N, 3)) + 1j * rng.standard_normal((N, 3))
FMM_ARGS = dict(fmm_maxdirect=10, device="none")
KERNELS = [("laplace", None), ("helmholtz", 3.0), ("helmholtz", 3.0 + 1.0j)]


def direct(tgt, src, q, kappa=0, normals=None, gradient=False, curl=False):
    """sum_i G(x-y_i) q_i, G = exp(i kappa r)/(4 pi r); pairs with r = 0 are skipped"""
    d = tgt[:, None, :] - src[None, :, :]
    r = np.linalg.norm(d, axis=2)
    ok = r > 0
    r = np.where(ok, r, 1.0)
    G = ok * np.exp(1j * kappa * r) / (4 * np.pi * r)
    dG = ((1j * kappa * r - 1) * G / r**2)[:, :, None] * d
    if curl:
        return np.cross(dG, q[None]).sum(axis=1)
    if normals is not None:
        return -np.einsum("jik,ik->ji", dG, normals) @ q
    return np.einsum("jik,i->jk", dG, q) if gradient else G @ q


def check(u, ref, tol=1e-6):
    assert u.shape == ref.shape
    assert np.max(np.abs(u - ref)) / np.max(np.abs(ref)) < tol


@pytest.mark.parametrize("kernel, kappa", KERNELS)
def test_charges(kernel, kappa):
    q = Q.real if kernel == "laplace" else Q
    ref = direct(TGT, SRC, q, kappa or 0)
    check(FMM(SRC, TGT, kernel=kernel, kappa=kappa, **FMM_ARGS) @ q, ref.real if kernel == "laplace" else ref)


def test_dipoles_gradient_and_coincident_targets():
    op = FMM(SRC, TGT, kernel="helmholtz", kappa=3.0, source="dipoles", source_normals=NRM, **FMM_ARGS)
    check(op @ Q, direct(TGT, SRC, Q, 3.0, normals=NRM))
    op = FMM(SRC, TGT, kernel="helmholtz", kappa=3.0, target="gradient", **FMM_ARGS)
    check(op @ Q, direct(TGT, SRC, Q, 3.0, gradient=True))
    op = FMM(SRC, kernel="helmholtz", kappa=3.0, **FMM_ARGS)
    check(op @ Q, direct(SRC, SRC, Q, 3.0))


def test_components_and_curl():
    kw = dict(kernel="helmholtz", kappa=3.0, **FMM_ARGS)
    u = FMM(SRC, TGT, components=3, **kw) @ Q3
    assert u.shape == (M, 3)
    check(u, np.stack([FMM(SRC, TGT, **kw) @ Q3[:, c] for c in range(3)], axis=1), 1e-12)
    check(FMM(SRC, TGT, source="curl", components=3, **kw) @ Q3, direct(TGT, SRC, Q3, 3.0, curl=True))


def test_fp32():
    q = Q.real
    check(FMM(SRC, TGT, fp32=True, **FMM_ARGS) @ q, direct(TGT, SRC, q).real, 1e-5)


def test_basematrix():
    op = FMM(SRC, TGT, target="gradient", **FMM_ARGS)
    assert op.shape == (3 * M, N)
    x = op.mat.CreateRowVector()
    x.FV().NumPy()[:] = Q.real
    y = op.mat.CreateColVector()
    y.data = op.mat * x
    check(y.FV().NumPy(), (op @ Q.real).ravel(), 1e-12)


def test_singular_expansion():
    mp = SingularMLExpansion((0, 0, 0), 1.0, 3.0)
    src, x = SRC - 0.5, rng.standard_normal((50, 3))
    x *= 3.5 / np.linalg.norm(x, axis=1)[:, None]
    mp.AddCharges(src, Q)
    mp.Calc()
    check(mp.Evaluate(x), direct(x, src, Q, 3.0))


def test_module_layout():
    assert bem.SingularMLExpansion is fmm.SingularMLExpansion
    assert not hasattr(bem, "FMM")
