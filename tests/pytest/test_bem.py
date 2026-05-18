
import pytest
import numpy as np

from netgen.occ import *
from ngsolve import *
from ngsolve.bem import *


@pytest.mark.parametrize(
    "order_source_1, order_source_2, order_target_1, order_target_2",
    [
        (80, 80, 80, 80),
        (80, 82, 79, 81),
    ],
)
def test_transform(order_source_1, order_source_2, order_target_1, order_target_2):
    box = Box((-10,-10,-10), (10,10,10))
    mesh = Mesh(OCCGeometry(box).GenerateMesh(maxh=5))

    kappa = 20
    S = SingularExpansionCF(order_source_1, kappa, (0,0,0), rad=1)
    S.AddCharge((0.3, -0.1,0.5), 1)

    dx,dy,dz = 0.1, 0.2, 0.3
    S2 = SingularExpansionCF(order_source_2, kappa, (dx,dy,dz), rad=2)
    S.Transform(S2)

    meshpnt = mesh(1,1,1)
    assert S(meshpnt) == pytest.approx(S2(meshpnt), rel=1e-10)

    R = RegularExpansionCF(order_target_1, kappa, (3,3,3), rad=2)
    S2.Transform(R)
    meshpnt = mesh(2,2,2)
    assert S(meshpnt) == pytest.approx(R(meshpnt), rel=1e-10)
    
    R2 = RegularExpansionCF(order_target_2, kappa, (3+dx,3+dy,3+dz), rad=1)
    R.Transform(R2)
    meshpnt = mesh(3,3,3)
    assert S(meshpnt) == pytest.approx(R2(meshpnt), rel=1e-10)


def test_singularmltransform():
    box = Box((-10,-10,-10), (10,10,10))
    mesh = Mesh(OCCGeometry(box).GenerateMesh(maxh=5))

    kappa = 0.01
    order = 20
    S = SingularMLExpansionCF((0,0,0), r=1, kappa=kappa)
    num = 200
    for i in range(num):
        z = i/num
        S.expansion.AddCharge((0.1, 0, z), 1/num)

    val1 = S(mesh(1,1,4))
    S.expansion.Calc()
    val2 = S(mesh(1,1,4))

    assert val1 == pytest.approx(val2)


def _sphere_mesh(maxh=0.25):
    sphere = Sphere((0, 0, 0), 1)
    return Mesh(OCCGeometry(sphere).GenerateMesh(maxh=maxh)).Curve(1)


def _fmm_test_mesh(geometry):
    if geometry == "sphere":
        return _sphere_mesh(maxh=0.25)
    if geometry == "box":
        box = Box((-1, -1, -1), (1, 1, 1))
        return Mesh(OCCGeometry(box).GenerateMesh(maxh=0.4)).Curve(1)
    if geometry == "quad_sphere":
        sp = Glue(Sphere((0, 0, 0), 1).faces)
        return Mesh(OCCGeometry(sp).GenerateMesh(maxh=0.25, quad_dominated=True)).Curve(3)
    raise ValueError(geometry)


def test_laplace_sl_calc_submatrix_matches_dense():
    mesh = _sphere_mesh(maxh=0.4)
    fes = SurfaceL2(mesh, order=1, complex=True)
    u, v = fes.TnT()

    op = LaplaceSL(u * ds) * v * ds

    rows = np.asarray([0, 5, 10, 15], dtype=np.int32)
    cols = np.asarray(
        [fes.ndof - 20, fes.ndof - 15, fes.ndof - 10, fes.ndof - 5],
        dtype=np.int32,
    )
    with TaskManager():
        mat = op.mat
        dense = np.asarray(mat.ToDense())

    submatrix = np.asarray(op.CalcSubMatrix(rows, cols))
    reference = dense[np.ix_(rows, cols)]

    assert submatrix.shape == (len(rows), len(cols))
    assert np.all(np.isfinite(submatrix))
    np.testing.assert_allclose(submatrix, reference, rtol=1e-8, atol=1e-10)


def _maxwell_sl(u, v, kappa, **kwargs):
    return kappa * HelmholtzSL(
        u.Trace() * ds, kappa, **kwargs
    ) * v.Trace() * ds - 1 / kappa * HelmholtzSL(
        div(u.Trace()) * ds, kappa, **kwargs
    ) * div(v.Trace()) * ds


def _fmm_direct_operators(operator_name, mesh):
    kappa = 1.5
    fmm_kwargs = {"use_fmm": True}
    higher_order_fmm_kwargs = {"use_fmm": True, "fmm_minorder": 30}

    if operator_name in (
        "LaplaceSL",
        "LaplaceDL",
        "HelmholtzSL",
        "HelmholtzDL",
        "HelmholtzCF",
    ):
        fes = SurfaceL2(mesh, order=1, complex=True)
        u, v = fes.TnT()
        if operator_name == "LaplaceSL":
            return (
                LaplaceSL(u * ds, **fmm_kwargs) * v * ds,
                LaplaceSL(u * ds, use_fmm=False) * v * ds,
            )
        if operator_name == "LaplaceDL":
            return (
                LaplaceDL(u * ds, **higher_order_fmm_kwargs) * v * ds,
                LaplaceDL(u * ds, use_fmm=False) * v * ds,
            )
        if operator_name == "HelmholtzSL":
            return (
                HelmholtzSL(u * ds, kappa, **fmm_kwargs) * v * ds,
                HelmholtzSL(u * ds, kappa, use_fmm=False) * v * ds,
            )
        if operator_name == "HelmholtzDL":
            return (
                HelmholtzDL(u * ds, kappa, **higher_order_fmm_kwargs) * v * ds,
                HelmholtzDL(u * ds, kappa, use_fmm=False) * v * ds,
            )
        return (
            HelmholtzCF(u * ds, kappa, **higher_order_fmm_kwargs) * v * ds,
            HelmholtzCF(u * ds, kappa, use_fmm=False) * v * ds,
        )

    if operator_name == "LameSL":
        fes = VectorH1(mesh, order=1)
        u, v = fes.TnT()
        return (
            LameSL(u * ds, E=2.0, nu=0.25, **fmm_kwargs) * v * ds,
            LameSL(u * ds, E=2.0, nu=0.25, use_fmm=False) * v * ds,
        )

    if operator_name == "MaxwellSL":
        hdiv = HDivSurface(mesh, order=1, complex=True)
        u, v = hdiv.TnT()
        return (
            _maxwell_sl(u, v, kappa, **fmm_kwargs),
            _maxwell_sl(u, v, kappa, use_fmm=False),
        )

    if operator_name == "MaxwellDL":
        hdiv = HDivSurface(mesh, order=1, complex=True)
        hcurl = HCurl(mesh, order=1, complex=True)
        u = hcurl.TrialFunction()
        v = hdiv.TestFunction()
        return (
            MaxwellDL(
                u.Operator("rotated_trace") * ds, kappa, **higher_order_fmm_kwargs
            ) * v.Trace() * ds,
            MaxwellDL(
                u.Operator("rotated_trace") * ds, kappa, use_fmm=False
            ) * v.Trace() * ds,
        )

    raise ValueError(operator_name)


@pytest.mark.parametrize("geometry", ["sphere", "box", "quad_sphere"])
@pytest.mark.parametrize(
    "operator_name",
    [
        "LaplaceSL",
        "LaplaceDL",
        "HelmholtzSL",
        "HelmholtzDL",
        "HelmholtzCF",
        "LameSL",
        "MaxwellSL",
        "MaxwellDL",
    ],
)
def test_fmm_and_direct_matrix_action(operator_name, geometry):
    mesh = _fmm_test_mesh(geometry)
    op_fmm, op_direct = _fmm_direct_operators(operator_name, mesh)

    with TaskManager():
        mat_fmm = op_fmm.mat
        mat_direct = op_direct.mat

    x = mat_fmm.CreateRowVector()
    x.FV().NumPy()[:] = np.linspace(0.25, 1.25, mat_fmm.width)
    y_fmm = mat_fmm.CreateColVector()
    y_direct = mat_direct.CreateColVector()
    y_fmm.data = mat_fmm * x
    y_direct.data = mat_direct * x

    diff = y_fmm.FV().NumPy() - y_direct.FV().NumPy()
    relerr = np.linalg.norm(diff) / np.linalg.norm(y_direct.FV().NumPy())
    assert relerr < 1e-8


@pytest.mark.parametrize("operator_name", ["LaplaceSL", "HelmholtzSL"])
def test_potential_operator_local_expansion_matches_direct_potential(operator_name):
    mesh = _sphere_mesh()
    fes = SurfaceL2(mesh, order=1, complex=operator_name == "HelmholtzSL")
    u = fes.TrialFunction()
    gfu = GridFunction(fes)
    gfu.vec.FV().NumPy()[:] = np.linspace(1, 2, fes.ndof)
    if operator_name == "LaplaceSL":
        potop = LaplaceSL(u * ds)
    else:
        potop = HelmholtzSL(u * ds, kappa=1.5)

    screen = WorkPlane(Axes((3, 0, 0), Y, Z)).RectangleC(0.5, 0.5).Face()
    screen.faces.name = "screen"
    target_mesh = Mesh(OCCGeometry(screen).GenerateMesh(maxh=0.4)).Curve(1)
    target_boundary = target_mesh.Boundaries("screen")

    direct = potop(gfu)
    local = potop(gfu, target_boundary)

    error = Integrate(local - direct, target_mesh, definedon=target_boundary)
    assert abs(error) < 1e-12


# Same setup as:
# https://github.com/Weggler/docu-ngsbem/blob/main/convergence_timing/Laplace_DtN_Convergence.py
def test_laplace_dtn_cg_solve_converges_to_docu_fundamental_solution():
    exact = 1 / sqrt((x - 1) ** 2 + (y - 1) ** 2 + (z - 1) ** 2)
    exact_grad = CF((exact.Diff(x), exact.Diff(y), exact.Diff(z)))
    exact_neumann = exact_grad * specialcf.normal(3)
    errors = []
    iterations = []

    for maxh in [0.45, 0.35, 0.28]:
        sphere = Sphere((0, 0, 0), 1)
        mesh = Mesh(OCCGeometry(sphere).GenerateMesh(maxh=maxh))
        order = 2
        fes_l2 = SurfaceL2(mesh, order=order - 1, dual_mapping=True)
        u, v = fes_l2.TnT()
        fes_h1 = H1(mesh, order=order)
        u_h1, v_h1 = fes_h1.TnT()
        dirichlet = GridFunction(fes_h1)
        dirichlet.Interpolate(exact)
        neumann = GridFunction(fes_l2)

        with TaskManager():
            pre = BilinearForm(u * v * ds, diagonal=True).Assemble().mat.Inverse()
            op = LaplaceSL(u * ds) * v * ds
            mass = BilinearForm(u_h1 * v * ds).Assemble()
            double_layer = LaplaceDL(u_h1 * ds) * v * ds
            rhs = ((0.5 * mass.mat + double_layer.mat) * dirichlet.vec).Evaluate()
            inv = solvers.CGSolver(op.mat, pre, tol=1e-10, maxiter=200, printrates=False)
            neumann.vec.data = inv * rhs

        error = sqrt(
            Integrate((exact_neumann - neumann) ** 2, mesh.Boundaries(".*"), BND)
        )
        errors.append(float(error))
        iterations.append(inv.iterations)

    assert errors[2] < errors[1] < errors[0]
    assert 0 < max(iterations) < 80
