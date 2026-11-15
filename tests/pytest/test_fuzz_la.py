import gc, pickle, subprocess, sys, textwrap
import numpy as np
import pytest
from ngsolve import *
from ngsolve.meshes import MakeStructured2DMesh


@pytest.fixture(scope="module")
def mesh():
    return MakeStructured2DMesh(nx=2, ny=2)


def run_script(code):
    return subprocess.run([sys.executable, "-c", textwrap.dedent(code)],
                          capture_output=True, text=True, timeout=120)


@pytest.mark.parametrize("inverse", ["sparsecholesky", "umfpack", "accelerate"])
def test_inverse_update_after_matrix_died(mesh, inverse):
    fes = H1(mesh, order=2)
    u, v = fes.TnT()
    a = BilinearForm(u*v*dx).Assemble()
    try:
        inv = a.mat.Inverse(fes.FreeDofs(), inverse=inverse)
    except Exception:
        pytest.skip(inverse + " not available")
    x = a.mat.CreateColVector(); x[:] = 1
    y0 = (inv * x).Evaluate()
    del a; gc.collect()
    with pytest.raises(Exception, match="does not exist anymore"):
        inv.Update()
    if inverse == "accelerate":
        return
    # the factorization is still usable
    assert np.allclose((inv * x).Evaluate().FV().NumPy(), y0.FV().NumPy())


def test_matrix_size_exceeds_int():
    with pytest.raises(ValueError):
        IdentityMatrix(2**31)
    r = la.Real2ComplexMatrix(IdentityMatrix(3))
    x = BaseVector(3, complex=True); x[:] = 1j
    assert np.allclose((r * x).Evaluate().FV().NumPy(), 1j)


def test_blockvector_empty():
    b = la.BlockVector([])
    assert len(b) == 0
    assert not b.is_complex
    try:
        np.asarray(b)
    except Exception:
        pass
    with pytest.raises(Exception):
        la.BlockVector([None])


def test_multivector_empty_orthogonalize():
    mv = MultiVector(BaseVector(5), 0)
    mv.Orthogonalize()
    R = mv.Orthogonalize(IdentityMatrix(5))
    assert R.h == 0 and R.w == 0
    mv = MultiVector(BaseVector(5), 2)
    mv[0][:] = 1; mv[1][:] = 0; mv[1][0] = 1
    mv.Orthogonalize()
    assert abs(InnerProduct(mv[0], mv[1])) < 1e-12


def test_eig_needs_square_matrix(mesh):
    v = FacetFESpace(mesh, order=1).TestFunction()
    with pytest.raises(Exception, match="square"):
        v.Operator("dual").Eig()


def test_proxy_dt_missing_operator(mesh):
    v = NumberSpace(mesh).TestFunction()
    assert v.Deriv() is None
    assert v.dt.Deriv() is None
    with pytest.raises(Exception):
        v.dt.Operator("Grad")
    w = HCurlCurl(mesh).TestFunction()
    assert list(Grad(w.dt).dims) == list(Grad(w).dims)
    u = H1(mesh, order=1).TrialFunction()
    assert list(u.dt.Deriv().dims) == list(u.Deriv().dims)


def test_frozen_missing_operator(mesh):
    v = NumberSpace(mesh).TestFunction()
    with pytest.raises(Exception):
        Grad(v.Freeze())
    u = H1(mesh, order=1).TrialFunction()
    assert list(Grad(u.Freeze()).dims) == list(Grad(u).dims)


@pytest.mark.parametrize("space", [HDivDiv, TangentialFacetFESpace])
def test_diffshape_proxy_direction(mesh, space):
    u, v = space(mesh, order=1).TnT()
    with pytest.raises(Exception):
        u.DiffShape(direction=v)
    w = VectorH1(mesh, order=1).TestFunction()
    u.DiffShape(direction=w)


def test_generate_l2_code_order():
    with pytest.raises(ValueError):
        fem.GenerateL2ElementCode(-2)
    assert "float" in fem.GenerateL2ElementCode(2)


def test_placeholder_pickle(mesh):
    p = fem.PlaceholderCF(x*y)
    p2 = pickle.loads(pickle.dumps(p))
    assert type(p2) is type(p)
    assert abs(Integrate(p2, mesh) - 0.25) < 1e-12


TAMPER = """
    import pickle, re
    from ngsolve import *
    b = bytearray(pickle.dumps(Conj(log(x))))
    # bytes inside class names are handled by netgen
    names = [range(m.start(), m.end()) for m in re.finditer(rb"ngfem::[A-Za-z_:<>]+", bytes(b))]
    for i in range(len(b)):
        if SKIP_NAMES and any(i in r for r in names): continue
        bb = bytearray(b); bb[i] = (bb[i] + 44) % 256
        try: pickle.loads(bytes(bb))
        except Exception: pass
    for n in range(len(b)):
        try: pickle.loads(bytes(b[:n]))
        except Exception: pass
    print("done")
    """


def test_tampered_pickle():
    r = run_script("SKIP_NAMES = True\n" + textwrap.dedent(TAMPER))
    assert r.returncode == 0 and "done" in r.stdout, r.stderr[-2000:]


@pytest.mark.xfail(reason="netgen Archive: unknown class name gives a null creator (archive.hpp)")
def test_tampered_pickle_class_name():
    r = run_script("SKIP_NAMES = False\n" + textwrap.dedent(TAMPER))
    assert r.returncode == 0 and "done" in r.stdout, r.stderr[-2000:]


@pytest.mark.parametrize("expr", [
    "BaseVector(10**12)",
    "BaseVector(2**40, complex=True)",
    "BaseVector(2**62, entrysize=8)",
    "MultiVector(BaseVector(10), 10**10)",
    "MultiVector(10**6, 10**8, False)",
    "BaseVector(5).CreateVectors(10**9)",
    "Matrix(10**6, 10**6)",
    "Matrix(2**31-1, 2**31-1, complex=True)",
])
def test_huge_sizes_memoryerror(expr):
    with pytest.raises(MemoryError):
        eval(expr)


@pytest.mark.parametrize("expr", [
    "la.SparseMatrixd.CreateFromCOO([0],[0],[1.0], 10**12, 10**12)",
    "IntegrationRule(ET.TET, 1000)",
    "BSpline(10**7, [0.5], [5])",
    "x.Compile(realcompile=True, maxderiv=-10**9)",
    "x.Compile(maxderiv=3)",
])
def test_huge_sizes_valueerror(expr):
    if "SparseMatrixd" in expr and not hasattr(la, "SparseMatrixd"):
        pytest.skip("SparseMatrix class name differs")
    with pytest.raises(ValueError):
        eval(expr)


def test_allocation_error_is_ngexception():
    from netgen.libngpy._meshing import NgException
    with pytest.raises(NgException):
        BaseVector(10**12)


def test_large_sizes_still_work():
    assert len(BaseVector(10**6)) == 10**6
    assert len(MultiVector(BaseVector(10), 1000)) == 1000
    assert len(BaseVector(5).CreateVectors(100)) == 100
    assert Matrix(1000, 1000).h == 1000
    assert len(IntegrationRule(ET.TET, 40)) > 0
    assert BSpline(5, [0, 1, 2, 3], [0, 1, 4, 9])(1.5) >= 0
