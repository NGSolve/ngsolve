import gc, pickle, subprocess, sys, textwrap
import numpy as np
import pytest
import scipy.sparse.linalg as spla
from ngsolve import *
from ngsolve.krylovspace import CG, GMRES, CGSolver, GMRESSolver
from ngsolve.eigenvalues import PINVIT, LOBPCG
from netgen.geom2d import unit_square


def diag_mat(n, val=2.):
    return la.SparseMatrixd.CreateFromCOO(list(range(n)), list(range(n)), [val]*n, n, n)


@pytest.fixture(scope="module")
def sys2d():
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.3))
    fes = H1(mesh, order=2, dirichlet=".*")
    u, v = fes.TnT()
    a = BilinearForm(grad(u)*grad(v)*dx + u*v*dx).Assemble()
    m = BilinearForm(u*v*dx).Assemble()
    return fes, a, m


def run_script(code):
    return subprocess.run([sys.executable, "-c", textwrap.dedent(code)],
                          capture_output=True, text=True, timeout=120)


def test_chebyshev_holds_operands():
    def f():
        return la.ChebyshevIteration(mat=diag_mat(1000), pre=IdentityMatrix(1000),
                                     steps=3, lam_min=0.5, lam_max=3)
    c = f(); gc.collect()
    junk = [BaseVector(1000) for i in range(50)]
    x = BaseVector(1000); y = BaseVector(1000); x[:] = 1
    y.data = c * x
    assert np.isfinite(y.Norm())


def test_asvector_keeps_matrix_alive():
    def f():
        A = diag_mat(500, 3.)
        return A.AsVector()
    vec = f(); gc.collect()
    junk = [BaseVector(vec.size) for i in range(5)]
    for j in junk: j[:] = 7
    vec[:] = 1
    assert vec.Norm() == pytest.approx(np.sqrt(500))


def test_register_inverse_bad_creator(sys2d):
    fes, a, m = sys2d
    la.RegisterInverseType("mybad_la", lambda mat: 17)
    with pytest.raises(Exception, match="SparseFactorizationInterface"):
        a.mat.Inverse(fes.FreeDofs(), inverse="mybad_la")


REG = """
from ngsolve import *
from netgen.geom2d import unit_square
class MyInv(la.SparseFactorizationInterface):
    def Analyze(self): pass
    def Factor(self): pass
    def Solve(self, b, x): x.data = b
la.RegisterInverseType("myinv", lambda m: MyInv(m))
"""


def test_register_inverse_exit_clean():
    r = run_script(REG)
    assert r.returncode == 0, r.stderr


def test_register_inverse_use_exit_clean():
    r = run_script(REG + """
mesh = Mesh(unit_square.GenerateMesh(maxh=0.3))
fes = H1(mesh, order=2, dirichlet=".*")
u,v = fes.TnT()
a = BilinearForm(grad(u)*grad(v)*dx + u*v*dx).Assemble()
inv = a.mat.Inverse(fes.FreeDofs(), inverse="myinv")
x = a.mat.CreateRowVector(); y = a.mat.CreateColVector(); x[:] = 1
y.data = inv * x
""")
    assert r.returncode == 0, r.stderr


def test_pylinearoperator_legit():
    op = spla.LinearOperator((5, 5), matvec=lambda x: 2*np.asarray(x), dtype=np.float64)
    m = BaseMatrix(op)
    x = BaseVector(5); y = BaseVector(5); x[:] = 1
    y.data = m * x
    assert y.Norm() == pytest.approx(2*np.sqrt(5))
    y.data += m * x
    assert y.Norm() == pytest.approx(4*np.sqrt(5))


def test_pylinearoperator_pylinop_in_thread_free_context():
    # Mult releases the GIL in the binding; the operator must reacquire it
    op = spla.LinearOperator((5, 5), matvec=lambda x: np.asarray(x), dtype=np.float64)
    m = BaseMatrix(op)
    x = BaseVector(5); y = BaseVector(5); x[:] = 1
    m.Mult(x, y)
    assert y.Norm() == pytest.approx(np.sqrt(5))


@pytest.mark.parametrize("n", [7, 700000, 2])
def test_pylinearoperator_bad_size(n):
    class Op:
        shape = (5, 5); dtype = np.dtype(np.float64)
        def __mul__(self, x): return np.ones(n)
    m = BaseMatrix(Op())
    x = BaseVector(5); y = BaseVector(5); x[:] = 1
    with pytest.raises(Exception, match="wrong size"):
        m.Mult(x, y)


def test_pylinearoperator_bad_type():
    class Op:
        shape = (5, 5); dtype = np.dtype(np.float64)
        def __mul__(self, x): return "nonsense"
    m = BaseMatrix(Op())
    x = BaseVector(5); y = BaseVector(5); x[:] = 1
    with pytest.raises(Exception):
        m.Mult(x, y)


def test_pickle_empty_multivector():
    x = BaseVector(10); x[:] = 1
    mv = MultiVector(x, 3)
    for i in range(3): mv[i][:] = i+1
    e = pickle.loads(pickle.dumps(MultiVector(x, 0)))
    assert len(e) == 0
    assert e.RefVec().size == 10 if hasattr(e, "RefVec") else True
    ip = e.InnerProduct(mv)
    assert np.asarray(ip).shape[0] == 3
    e.Extend(2)
    assert len(e) == 2
    mv2 = pickle.loads(pickle.dumps(mv))
    assert len(mv2) == 3 and mv2[1].Norm() == pytest.approx(mv[1].Norm())


def test_multivector_innerproduct_mismatch():
    mv_r = MultiVector(BaseVector(100), 2)
    mv_c = MultiVector(BaseVector(100, complex=True), 2)
    mv_small = MultiVector(BaseVector(10), 2)
    with pytest.raises(Exception):
        mv_c.InnerProduct(mv_r)
    with pytest.raises(Exception):
        mv_r.InnerProduct(mv_c)
    with pytest.raises(Exception):
        mv_r.InnerProduct(mv_small)
    with pytest.raises(Exception):
        mv_small.InnerProduct(BaseVector(100))
    assert np.asarray(mv_r.InnerProduct(mv_r)).shape == (2, 2)


def test_cg_nan_breakdown(sys2d):
    fes, a, m = sys2d
    pre = a.mat.CreateSmoother(fes.FreeDofs())
    s = CGSolver(mat=a.mat, pre=pre, maxiter=100000000, tol=1e-30, printrates=False)
    x = a.mat.CreateRowVector(); x[:] = float("nan")
    with pytest.raises(ArithmeticError, match="breakdown"):
        s.Solve(x)


def test_solver_argument_validation(sys2d):
    fes, a, m = sys2d
    pre = a.mat.CreateSmoother(fes.FreeDofs())
    x = a.mat.CreateRowVector(); x[:] = 1
    for kw in (dict(maxiter=0), dict(maxiter=-5), dict(tol=-1), dict(tol=0), dict(tol=float("nan")),
               dict(maxsteps=0), dict(precision=-1)):
        with pytest.raises(ValueError):
            CGSolver(mat=a.mat, pre=pre, printrates=False, **kw)
    with pytest.raises(ValueError):
        CGSolver(mat=a.mat, printrates=False)
    for r in (0, -1):
        with pytest.raises(ValueError):
            GMRESSolver(mat=a.mat, pre=pre, restart=r, printrates=False)
    with pytest.raises(ValueError):
        GMRES(a.mat, x, pre=pre, restart=0, printrates=False)
    with pytest.raises(ValueError):
        CG(a.mat, x, pre=pre, maxsteps=0, printrates=False)
    # valid use still works
    y = CG(a.mat, x, pre=pre, maxsteps=50, printrates=False)
    free = np.array(fes.FreeDofs())
    r = x.CreateVector(); r.data = a.mat * y - x
    res = r.FV().NumPy()[free]
    assert np.linalg.norm(res) < 1e-6 * np.linalg.norm(x.FV().NumPy()[free])


def test_eigen_num_validation(sys2d):
    fes, a, m = sys2d
    pre = a.mat.CreateSmoother(fes.FreeDofs())
    short = BitArray(3); short.Set()
    with pytest.raises(ValueError):
        PINVIT(a.mat, m.mat, pre, num=500, maxit=2, printrates=False)
    with pytest.raises(ValueError):
        LOBPCG(a.mat, m.mat, pre, num=500, maxit=2, printrates=False)
    with pytest.raises(ValueError):
        PINVIT(a.mat, m.mat, Projector(short, False), num=3, maxit=3, printrates=False)
    with pytest.raises(ValueError):
        PINVIT(a.mat, m.mat, pre, num=0, printrates=False)
    lams, vecs = PINVIT(a.mat, m.mat, pre, num=2, maxit=20, printrates=False)
    assert len(lams) == 2
