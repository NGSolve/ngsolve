
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


@pytest.mark.parametrize("operator_name", ["LaplaceSL", "HelmholtzSL"])
def test_volume_sl_direct_and_fmm_match(operator_name):
    mesh = Mesh(unit_cube.GenerateMesh(maxh=2))
    is_helmholtz = operator_name == "HelmholtzSL"
    fes = L2(mesh, order=0, complex=is_helmholtz)
    u, v = fes.TnT()
    dxi = dx(bonus_intorder=2)

    if is_helmholtz:
        direct = HelmholtzSL(u * dxi, 1.2 + 0.3j, use_fmm=False) * v * dxi
        fast = HelmholtzSL(u * dxi, 1.2 + 0.3j, use_fmm=True) * v * dxi
    else:
        direct = LaplaceSL(u * dxi, use_fmm=False) * v * dxi
        fast = LaplaceSL(u * dxi, use_fmm=True) * v * dxi

    mat_direct = np.asarray(direct.mat.ToDense())
    mat_fast = np.asarray(fast.mat.ToDense())

    assert np.all(np.isfinite(mat_direct))
    np.testing.assert_allclose(mat_fast, mat_direct, rtol=1e-8, atol=1e-11)

    if is_helmholtz:
        ids = np.asarray([0, 3, 7], dtype=np.int32)
        submatrix = np.asarray(direct.CalcSubMatrix(ids, ids))
        np.testing.assert_allclose(
            submatrix, mat_direct[np.ix_(ids, ids)], rtol=1e-10, atol=1e-12
        )


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


def _fmm_direct_operators(operator_name, mesh, kappa=1.5, order=20, fmm_options=None):
    fmm_kwargs = {"use_fmm": True, "fmm_minorder": order}
    if fmm_options:
        fmm_kwargs.update(fmm_options)

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
                LaplaceDL(u * ds, **fmm_kwargs) * v * ds,
                LaplaceDL(u * ds, use_fmm=False) * v * ds,
            )
        if operator_name == "HelmholtzSL":
            return (
                HelmholtzSL(u * ds, kappa, **fmm_kwargs) * v * ds,
                HelmholtzSL(u * ds, kappa, use_fmm=False) * v * ds,
            )
        if operator_name == "HelmholtzDL":
            return (
                HelmholtzDL(u * ds, kappa, **fmm_kwargs) * v * ds,
                HelmholtzDL(u * ds, kappa, use_fmm=False) * v * ds,
            )
        return (
            HelmholtzCF(u * ds, kappa, **fmm_kwargs) * v * ds,
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
                u.Operator("rotated_trace") * ds, kappa, **fmm_kwargs
            ) * v.Trace() * ds,
            MaxwellDL(
                u.Operator("rotated_trace") * ds, kappa, use_fmm=False
            ) * v.Trace() * ds,
        )

    raise ValueError(operator_name)

_fmm_operator_cases = [
    ("LaplaceSL", 1.5),
    ("LaplaceDL", 1.5),
    ("HelmholtzSL", 1.5),
    ("HelmholtzSL", 7.5 + 10j),
    ("HelmholtzDL", 1.5),
    ("HelmholtzDL", 7.5 + 10j),
    ("HelmholtzCF", 1.5),
    ("HelmholtzCF", 1.0 + 5.0j),
    ("LameSL", 1.5),
    ("MaxwellSL", 1.5),
    ("MaxwellDL", 1.5),
]

_fmm_matrix_action_cases = [
    (operator_name, kappa, 30, geometry)
    for operator_name, kappa in _fmm_operator_cases
    for geometry in ["sphere", "box", "quad_sphere"]
] + [
    (operator_name, kappa, 50, "sphere")
    for operator_name, kappa in [
        ("HelmholtzSL", 7.5 + 10j),
        ("HelmholtzDL", 7.5 + 10j),
        ("HelmholtzCF", 1.0 + 5.0j),
    ]
]


@pytest.mark.parametrize("operator_name, kappa, order, geometry", _fmm_matrix_action_cases)
def test_fmm_and_direct_matrix_action(operator_name, geometry, kappa, order):
    mesh = _fmm_test_mesh(geometry)
    op_fmm, op_direct = _fmm_direct_operators(operator_name, mesh, kappa, order)

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
    tol = 1e-7 if operator_name == "MaxwellDL" and geometry == "quad_sphere" else 1e-8
    assert relerr < tol


@pytest.mark.parametrize(
    "operator_name, kappa, order, maxdirect",
    [
        ("LaplaceSL", 1.5, 20, 5),
        ("LaplaceDL", 1.5, 20, 5),
        ("HelmholtzSL", 1.5, 20, 5),
        ("HelmholtzDL", 1.5, 20, 5),
        ("HelmholtzCF", 8.0, 20, 5),
        ("MaxwellDL", 1.5, 20, 5),
        ("LameSL", 1.5, 20, 5),
        ("HelmholtzCF", 100.0, 20, 25),
    ] + [
        (operator_name, kappa, 30, 25)
        for operator_name, kappa in [
            ("HelmholtzSL", 1.5),
            ("HelmholtzDL", 1.5),
            ("MaxwellDL", 1.5),
            ("HelmholtzSL", 7.5 + 10j),
            ("HelmholtzDL", 7.5 + 10j),
            ("HelmholtzCF", 1.0 + 5.0j),
        ]
    ],
)
def test_fp32_fmm_matches_fp64_far_field(operator_name, kappa, order, maxdirect):
    mesh = _sphere_mesh(maxh=0.4)

    with TaskManager():
        op64, _ = _fmm_direct_operators(
            operator_name, mesh, kappa, order=order,
            fmm_options={"fmm_maxdirect": maxdirect},
        )
        op32, _ = _fmm_direct_operators(
            operator_name, mesh, kappa, order=order,
            fmm_options={"fmm_maxdirect": maxdirect, "fp32": True},
        )

    info64 = op64.GetFMMInfo()
    info32 = op32.GetFMMInfo()
    assert info32["total_multipole_coefficients"] == info64["total_multipole_coefficients"]
    assert info32["multipole_memory_mb"] == pytest.approx(0.5 * info64["multipole_memory_mb"])
    assert info32["num_s2r"] > 0
    assert info32["direct_fallback_fraction"] < 1.0

    vec = op64.mat.CreateRowVector()
    expected_dtype = np.float64 if operator_name == "LameSL" else np.complex128
    assert vec.FV().NumPy().dtype == expected_dtype
    values = np.linspace(0.25, 1.25, op64.mat.width)
    vec.FV().NumPy()[:] = values
    if np.iscomplexobj(vec.FV().NumPy()):
        vec.FV().NumPy()[:] += 0.2j * values[::-1]
    out64 = op64.mat.CreateColVector()
    out32 = op32.mat.CreateColVector()
    with TaskManager():
        out64.data = op64.mat * vec
        out32.data = op32.mat * vec

    diff = out32.FV().NumPy() - out64.FV().NumPy()
    assert np.all(np.isfinite(out32.FV().NumPy()))
    relerr = np.linalg.norm(diff) / np.linalg.norm(out64.FV().NumPy())
    assert relerr < 1e-5

    reference32 = out32.FV().NumPy().copy()
    for _ in range(3):
        repeated = op32.mat.CreateColVector()
        with TaskManager():
            repeated.data = op32.mat * vec
        repeat_diff = repeated.FV().NumPy() - reference32
        repeat_relerr = np.linalg.norm(repeat_diff) / np.linalg.norm(reference32)
        assert repeat_relerr < 1e-5


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
