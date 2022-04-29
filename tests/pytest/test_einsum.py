#!/usr/bin/env python
# coding: utf-8

from netgen.csg import unit_cube
from ngsolve import *
import numpy as np
import pytest

mesh = Mesh(unit_cube.GenerateMesh())

fes = VectorH1(mesh)
u = GridFunction(fes)
F = Grad(u) + Id(3)
u.Interpolate((x ** 2 * z * y, y ** 3 * x, z ** 2 * x * y))

X0 = mesh(0.1, 0.1, 0.1)

p_ = [Parameter(float((F * F)(X0)[i])) for i in range(9)]
pF = CoefficientFunction(tuple(p_), dims=(3, 3)).MakeVariable()


def DetES(A, **options):
    return fem.Einsum('ijk,i,j,k->', fem.LeviCivitaSymbol(3), *[A[i, :] for i in range(3)], **options)


def TraceES(A, **options):
    return fem.Einsum('ii->', A, **options)


def TransposeES(A, **options):
    return fem.Einsum('ij->ji', A, **options)


def InvES(A, **options):
    _Det = DetES(A, **options)
    return TransposeES(_Det.Diff(A.MakeVariable()) / _Det, **options)


def op_counts(acf):
    steps = str(acf).split("\n")
    ops = {}
    for step in steps[1:]:
        if len(step) == 0 or step[0] == " ":
            continue
        op = str(step[step.find(":") + 2:])
        if op not in ops.keys():
            ops[op] = 0
        ops[op] += 1
    return ops


def same(cf1, cf2, tol=1e-12):
    if np.array(cf1(X0)).size != np.array(cf2(X0)).size:
        return False
    if np.max(np.abs(np.array(cf1(X0)) - np.array(cf2(X0)))) < tol:
        return True
    else:
        print(np.max(np.abs(np.array(cf1(X0)) - np.array(cf2(X0)))), " >= ", tol)
        return False


@pytest.mark.parametrize("use_legacy_ops", (True, False))
def test_blas(use_legacy_ops):
    options = {"use_legacy_ops": use_legacy_ops, "optimize_path": False}
    
    def check_optimization(cf, options, legacy_str_lines):
        if options["use_legacy_ops"]:
            cflines = str(cf).splitlines()
            for k, v in legacy_str_lines.items():
                if cflines[k].count(v) != 1:
                    return False
        return True

    op = fem.Einsum('ij,ik->jk', pF, pF, **options)
    op_blas = pF.trans * pF
    assert same(op, op_blas)
    assert check_optimization(op, options, {0: "optimized node matrix-matrix multiply",
                                            1: "Matrix transpose"})

    op = fem.Einsum('ij,kj->ik', pF, pF, **options)
    op_blas = pF * pF.trans
    assert same(op, op_blas)
    assert check_optimization(op, options, {0: "optimized node matrix-matrix multiply",
                                            12: "Matrix transpose"})

    op = fem.Einsum('ij,jk->ik', pF, pF, **options)
    op_blas = pF * pF
    assert same(op, op_blas)
    assert check_optimization(op, options, {0: "optimized node matrix-matrix multiply"})

    op = fem.Einsum('ij,j->i', pF, pF[:, 0], **options)
    op_blas = pF * pF[:, 0]
    assert same(op, op_blas)
    assert check_optimization(op, options, {0: "optimized node matrix-vector multiply"})

    op = fem.Einsum('ij,i->j', pF, pF[:, 0], **options)
    op_blas = pF.trans * pF[:, 0]
    assert same(op, op_blas)
    assert check_optimization(op, options, {0: "optimized node matrix-vector multiply",
                                            1: "Matrix transpose"})

    op = fem.Einsum('ii', pF, **options)
    op_blas = Trace(pF)
    assert same(op, op_blas)
    assert check_optimization(op, options, {0: "optimized node trace"})

    op = fem.Einsum('ij->ji', pF, **options)
    op_blas = pF.trans
    assert same(op, op_blas)
    assert check_optimization(op, options, {0: "optimized node Matrix transpose"})

    pFpF = OuterProduct(pF, pF).Reshape((3, 3, 3, 3))
    op = fem.Einsum('ijkl->jilk', pFpF, **options)
    op_blas = pFpF.TensorTranspose((1, 0, 3, 2))
    assert same(op, op_blas)
    assert check_optimization(op, options, {0: "with optimized node tensor-transpose [ 1, 0, 3, 2 ]"})


def test_identity_optimizations():
    def check_optimization(cf, legacy_str_lines):
        cflines = str(cf).splitlines()
        for k, v in legacy_str_lines.items():
            if cflines[k].count(v) != 1:
                return False
        return True

    options = {"optimize_identities": True, "optimize_path": True}
    I = Id(3)
    II = Id([3, 3])

    op = fem.Einsum('ik,kj->ij', I, pF, **options)
    op_noopt = fem.Einsum('ik,kj->ij', I, pF, optimize_identities=False)
    op_opt = pF
    assert same(op, op_opt)
    assert same(op, op_noopt)
    assert check_optimization(op, {0: "EinsumCF ik,kj->ij with optimized node EinsumCF ij->ij"})

    op = fem.Einsum('ijkl,kl->ij', II, pF, **options)
    op_noopt = fem.Einsum('ijkl,kl->ij', II, pF, optimize_identities=False, sparse_evaluation=False)
    op_opt = pF
    assert same(op, op_opt)
    assert same(op, op_noopt)
    assert check_optimization(op, {0: "EinsumCF ijkl,kl->ij with optimized node EinsumCF ij->ij"})

    op = fem.Einsum('ii,kj->kj', I, pF, **options)
    op_opt = 3.0 * pF
    assert same(op, op_opt)
    assert check_optimization(op, {0: "EinsumCF ii,kj->kj with optimized node EinsumCF kj,i->kj"})

    op = fem.Einsum('ijkl,jl->ik', II, pF, **options)
    op_opt = Trace(pF) * I
    assert same(op, op_opt)
    assert check_optimization(op, {0: "EinsumCF ijkl,jl->ik with optimized node EinsumCF ll,ik->ik"})

    op = fem.Einsum('ii,ij->ij', Id(3), pF, **options)
    op_opt = pF
    assert same(op, op_opt)

    op = fem.Einsum('ii->ii', pF, **options)
    op_opt = pF
    assert same(op, op_opt)
    assert check_optimization(op, {0: "EinsumCF ii->ii with optimized node unary operation ' '"})


def test_expansion():
    options = {"expand_einsum": True}
    op1 = fem.Einsum('ik,kj->ij', 1 * pF, 2 * pF, **options)
    op2 = fem.Einsum('ij,jl->il', 3 * pF, op1, **options)
    op2_e = fem.Einsum('ij,jk,kl->il', 3 * pF, 1 * pF, 2 * pF, **options)
    assert same(op2, op2_e)
    op3 = fem.Einsum('ik,jl->ijkl', op1, op2, **options)
    op3_e = fem.Einsum('io,ok,jm,mn,nl->ijkl', 1 * pF, 2 * pF, 3 * pF, 1 * pF, 2 * pF, **options)
    assert same(op3, op3_e)
    op4 = fem.Einsum('ijkl,jl->ik', op3, op1, **options)
    op4_e = fem.Einsum('io,ok,jm,mn,nl,jp,pl->ik', 1 * pF, 2 * pF, 3 * pF, 1 * pF, 2 * pF, 1 * pF, 2 * pF, **options)
    assert same(op4, op4_e)


def test_path_optimization():
    options = {"expand_einsum": True,
               "optimize_path": True,
               "optimize_identities": False}
    op1 = fem.Einsum('ik,kj->ij', 1 * pF, 2 * pF)
    opt1 = fem.Einsum('ik,kj->ij', 1 * pF, 2 * pF, **options)
    assert same(op1, opt1)
    op2 = fem.Einsum('ij,jl->il', 3 * pF, op1)
    op3 = fem.Einsum('ik,jl->ijkl', op1, op2)
    opt3 = fem.Einsum('ik,jl->ijkl', op1, op2, **options)
    assert same(op3, opt3)
    op4 = fem.Einsum('ik,jl->ijkl', Id(3), op2)
    opt4 = fem.Einsum('ik,jl->ijkl', Id(3), op2, **options)
    assert same(op4, opt4)
    op5 = fem.Einsum('ijkl,jmln->ikmn', op4, op3)
    opt5 = fem.Einsum('ijkl,jmln->ikmn', op4, op3, **options)
    options["optimize_identities"] = True
    opt5b = fem.Einsum('ijkl,jmln->ikmn', op4, op3, **options)
    assert same(op5, opt5)
    assert same(op5, opt5b)
    assert same(op5, opt5.Compile())
    assert same(op5, opt5.Compile(realcompile=True, wait=True))


@pytest.mark.parametrize("options", ({"expand_einsum": True, "optimize_path": True, "optimize_identities": True},))
def test_diff(options):
    Cv = fem.Einsum('ki,kj->ij', pF, pF).MakeVariable()
    b = fem.Einsum('ij,jk,kl->il', F, InvES(Cv, optimize_path=False), TransposeES(F, optimize_path=False))
    Psi = TraceES(b, optimize_path=False) + log(DetES(b, optimize_path=False))
    PsiOpt = TraceES(b, **options) + log(DetES(b, **options))

    _Cv = (pF.trans * pF).MakeVariable()
    _b = F * Inv(_Cv) * F.trans
    Psi_legacy = Trace(_b) + log(Det(_b))

    assert same(DetES(Cv), Det(_Cv))
    assert same(InvES(Cv), Inv(_Cv))

    assert same(b, _b)
    assert same(TraceES(b), Trace(_b))
    assert same(DetES(b), Det(_b))

    assert same(Psi, Psi_legacy)

    DPsi = Psi.Diff(Cv)
    DDPsi = DPsi.Diff(Cv)

    DPsi_legacy = Psi_legacy.Diff(_Cv)
    DDPsi_legacy = DPsi_legacy.Diff(_Cv)

    assert same(DPsi, DPsi_legacy)
    assert same(DDPsi, DDPsi_legacy)
    assert same(DDPsi, PsiOpt.Diff(Cv).Diff(Cv).Compile())


def test_tensor_diff():
    options = {"optimize_path": True}

    Cv = fem.Einsum('ki,kj->ij', pF, pF).MakeVariable()
    b = fem.Einsum('ij,jk,kl->il', F, InvES(Cv), TransposeES(F, **options))
    assert same(Trace(Cv).Diff(Cv), fem.Einsum('ii', Cv, **options).Diff(Cv))
    assert same(Trace(b).Diff(Cv), fem.Einsum('ii', b, **options).Diff(Cv))

    assert same(Trace(pF).Diff(pF), TraceES(pF, **options).Diff(pF))
    assert same(Trace(b).Diff(Cv), TraceES(b, **options).Diff(Cv))

    assert same(Det(pF).Diff(pF), DetES(pF, **options).Diff(pF))
    assert same(Det(b).Diff(Cv), DetES(b, **options).Diff(Cv))

    assert same(Inv(pF).Diff(pF), InvES(pF, **options).Diff(pF))
    assert same(Inv(b).Diff(Cv), InvES(b, **options).Diff(Cv))

    assert same(pF.trans.Diff(pF), TransposeES(pF, **options).Diff(pF))
    assert same(b.trans.Diff(Cv), TransposeES(b, **options).Diff(Cv))

    A_extend = np.zeros((3, 3, 2, 2))
    A_extend[1:, :-1, :, :] = np.einsum('ik,jl->ijkl', np.eye(2), np.eye(2))
    cf_extend = CoefficientFunction(tuple(A_extend.flatten().tolist()), dims=A_extend.shape)
    pF_extend = pF[:2, :2].ExtendDimension((3, 3), (1, 0))
    pF_extend_ES = fem.Einsum('ijkl,kl->ij', cf_extend, pF[:2, :2])
    assert same(pF_extend.Diff(pF), pF_extend_ES.Diff(pF))

    b_extend = b[:2, :2].ExtendDimension((3, 3), (1, 0))
    b_extend_ES = fem.Einsum('ijkl,kl->ij', cf_extend, b[:2, :2])
    assert same(b_extend.Diff(Cv), b_extend_ES.Diff(Cv))

    AA = fem.Einsum('ik,jl->ijkl', pF, b)
    assert same(AA.TensorTranspose((2, 0, 3, 1)).Diff(pF), fem.Einsum("ijkl->kilj", AA).Diff(pF))


def test_ellipses():
    options = {"optimize_path": False}

    def check_signature(cf, expected_signature):
        return str(cf).splitlines()[0].count(" " + expected_signature + ",") == 1

    Cv = fem.Einsum('ki,kj->ij', pF, pF).MakeVariable()
    b = fem.Einsum('ij,jk,kl->il', F, InvES(Cv, **options), TransposeES(F, **options))
    bb = fem.Einsum('ij,kl->ijkl', b, b)
    cb = fem.Einsum('i,jk->ikj', CF(2).Reshape((1, )), b)
    assert check_signature(fem.Einsum('i...,j...k,k...->ijk', cb, bb, Cv, **options), 'iAB,jABk,kA->ijk')
    assert check_signature(fem.Einsum('i...,j...k,k...->ij...k', cb, bb, Cv, **options), 'iAB,jABk,kA->ijABk')
    assert check_signature(fem.Einsum('i,...kj,kj->i...k', CF(2).Reshape((1,)), bb, Cv, **options), 'i,ABkj,kj->iABk')


def test_zero_detection():
    def check_optimization(cf):
        return str(cf).splitlines()[0].count("ZeroCF") == 1

    AA = fem.Einsum('ik,jl->ijkl', pF, 0 * pF)
    assert same(AA, 0 * fem.Einsum('ik,jl->ijkl', pF, pF))
    assert check_optimization(AA)

    # DiffJacobi does its own "Zero optimization"
    AA = fem.Einsum('ik,jl->ijkl', pF, pF).Diff((2*pF).MakeVariable())
    assert np.max(np.abs(np.array(AA.dims) - [3] * 6)) < 1e-12
    assert str(AA) == "ZeroCoefficientFunction"


if __name__ == "__main__":
    test_blas(False)
    test_blas(True)
    test_identity_optimizations()
    test_expansion()
    test_path_optimization()
    test_diff({"expand_einsum": True,
               "optimize_path": True,
               "optimize_identities": True})
    test_tensor_diff()
    test_ellipses()
    test_zero_detection()
