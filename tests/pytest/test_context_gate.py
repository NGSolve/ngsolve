"""
A CoefficientFunction bound to a context (a form on a space, an integral over
a VorB of a mesh, a point evaluation) is checked once, before any element is
touched: proxies must fit the form's spaces, geometric CFs the mesh, region-wise
CFs the regions, complex values the scalar type.
"""
import pytest
import netgen.geom2d as g2d
from netgen.geom2d import unit_square
from ngsolve import *
from ngsolve.comp import Hidden, Reorder
from ngsolve.meshes import MakeStructured2DMesh, MakeStructured3DMesh


@pytest.fixture(scope="module")
def m():
    return Mesh(unit_square.GenerateMesh(maxh=0.3))


@pytest.fixture(scope="module")
def m2mat():
    geo = g2d.SplineGeometry()
    geo.AddRectangle((0, 0), (1, 1), leftdomain=1)
    geo.AddRectangle((0.3, 0.3), (0.6, 0.6), leftdomain=2, rightdomain=1)
    geo.SetMaterial(1, "a")
    geo.SetMaterial(2, "b")
    return Mesh(geo.GenerateMesh(maxh=0.3))


def nrm(a):
    return Norm(a.mat.AsVector())


def raises(func, match=None):
    with pytest.raises(Exception, match=match):
        func()


# ---------------------------------------------------------------- positive cases

def test_lowordering_space(m):
    V = H1(m, order=3)
    u, v = V.TnT()
    a = BilinearForm(V.lospace)
    a += grad(u) * grad(v) * dx
    a.Assemble()
    assert nrm(a) > 0 and a.mat.height == V.lospace.ndof


@pytest.mark.parametrize("wrap", ["periodic", "compress", "discontinuous", "hidden", "reorder"])
def test_wrapper_with_base_proxies(wrap):
    mesh = MakeStructured2DMesh(quads=False, nx=4, ny=4, periodic_x=True)
    V = H1(mesh, order=2)
    W = {"periodic": lambda: Periodic(V),
         "compress": lambda: Compress(V),
         "discontinuous": lambda: Discontinuous(V),
         "hidden": lambda: Hidden(V),
         "reorder": lambda: Reorder(V)}[wrap]()
    u, v = V.TnT()
    a = BilinearForm(W)
    a += (grad(u) * grad(v) + u * v) * dx
    a.Assemble()
    f = LinearForm(W)
    f += v * dx
    f.Assemble()
    assert a.mat.height == W.ndof and f.vec.size == W.ndof


def test_compound_components(m):
    X = H1(m, order=2) * HDiv(m, order=1) * NumberSpace(m)
    (u, s, lam), (v, t, mu) = X.TnT()
    a = BilinearForm(X)
    a += (grad(u) * grad(v) + s * t + div(s) * v + u * div(t) + lam * v + u * mu) * dx
    a.Assemble()
    assert nrm(a) > 0


def test_same_family_other_space(m):
    # proxies of a separately built H1 act on the elements of the form's space
    Vh = H1(m, order=4)
    u, v = Vh.TnT()
    V1 = H1(m, order=1)
    a = BilinearForm(V1)
    a += grad(u) * grad(v) * dx
    a.Assemble()
    assert a.mat.height == V1.ndof


def test_trace_on_bnd(m):
    D = HDiv(m, order=1)
    s, t = D.TnT()
    a = BilinearForm(D)
    a += s.Trace() * t.Trace() * ds
    a.Assemble()
    H = H1(m, order=2)
    u, v = H.TnT()
    b = BilinearForm(H)
    b += u * v * ds
    b.Assemble()
    assert nrm(a) > 0 and nrm(b) > 0


def test_element_boundary_and_skeleton(m):
    n = specialcf.normal(2)
    D = HDiv(m, order=1)
    s, t = D.TnT()
    a = BilinearForm(D)
    a += (s * n) * (t * n) * dx(element_boundary=True)
    a.Assemble()
    L = L2(m, order=1, dgjumps=True)
    u, v = L.TnT()
    b = BilinearForm(L)
    b += (u - u.Other()) * (v - v.Other()) * dx(skeleton=True)
    b += (grad(u) * n) * v * dx(element_boundary=True)
    b += u * v * ds(skeleton=True)
    b.Assemble()
    assert nrm(a) > 0 and nrm(b) > 0


def test_skeleton_without_dgjumps_apply_only(m):
    # explicit DG: matrix-free application does not need dgjumps
    L = L2(m, order=1)
    u, v = L.TnT()
    b = CF((1, 0.3))
    bn = b * specialcf.normal(2)
    a = BilinearForm(L)
    a += bn * IfPos(bn, u, u.Other()) * (v - v.Other()) * dx(skeleton=True)
    gf = GridFunction(L)
    gf.Set(x)
    r = gf.vec.CreateVector()
    a.Apply(gf.vec, r)
    assert Norm(r) > 0


def test_residual_term_in_bilinearform(m):
    # test-function only terms are residual contributions of nonlinear forms
    V = H1(m, order=2, dirichlet=".*")
    u, v = V.TnT()
    a = BilinearForm(V)
    a += (grad(u) * grad(v) + u**3 * v) * dx - 1 * v * dx
    gf = GridFunction(V)
    r = gf.vec.CreateVector()
    a.Apply(gf.vec, r)
    a.AssembleLinearization(gf.vec)
    assert Norm(r) > 0


def test_mixed_dgjumps_on_one_side(m):
    H = H1(m, order=2)
    L = L2(m, order=1, dgjumps=True)
    a = BilinearForm(trialspace=H, testspace=L)
    a += H.TrialFunction() * L.TestFunction() * dx
    a.Assemble()
    assert a.mat.height == L.ndof and a.mat.width == H.ndof


def test_normal_on_boundaries(m):
    n = specialcf.normal(2)
    assert Integrate(n * n, m, BND) == pytest.approx(4)
    m3 = MakeStructured3DMesh(nx=2, ny=2, nz=2)
    assert Integrate(specialcf.normal(3) * specialcf.normal(3), m3, BND) == pytest.approx(6)


def test_domainwise(m2mat):
    assert Integrate(CF([1, 2]), m2mat) == pytest.approx(1 + 0.09)
    assert Integrate(m2mat.MaterialCF({"b": 1}, default=0), m2mat) == pytest.approx(0.09)
    assert Integrate(CF([1, 2, 3]), m2mat) == pytest.approx(1 + 0.09)


def test_complex_forms(m):
    C = H1(m, order=2, complex=True)
    u, v = C.TnT()
    a = BilinearForm(C)
    a += (grad(u) * grad(v) + 1j * u * v) * dx
    a.Assemble()
    f = LinearForm(C)
    f += 1j * v * dx
    f.Assemble()
    R = H1(m, order=2)
    ur, vr = R.TnT()
    b = BilinearForm(R)
    b += (ur * vr * (1j * x).imag + Norm(1j * x) * ur * vr) * dx
    b.Assemble()
    assert nrm(a) > 0 and Norm(f.vec) > 0 and nrm(b) > 0


def test_interpolate_operator(m):
    H = H1(m, order=2)
    u, v = H.TnT()
    a = BilinearForm(H)
    a += Interpolate(grad(u), VectorL2(m, order=1)) * grad(v) * dx
    a.Assemble()
    # proxies of a separately built space of the same family
    s = H1(m, order=3).TrialFunction()
    b = BilinearForm(H)
    b += Interpolate(s, L2(m, order=0)) * v * dx
    b.Assemble()
    assert nrm(a) > 0 and nrm(b) > 0
    t = HDiv(m, order=1).TrialFunction()
    c = BilinearForm(H)
    raises(lambda: c.__iadd__(Interpolate(t, VectorL2(m, order=0)) * CF((1, 1)) * v * dx))


def test_point_evaluation(m):
    H = H1(m, order=2)
    gf = GridFunction(H)
    gf.Set(x * y)
    assert gf(m(0.3, 0.4)) == pytest.approx(0.12)
    cf = sin(x) * gf + y
    for i in range(3):
        assert cf(m(0.3, 0.4)) == pytest.approx(sin(0.3) * 0.12 + 0.4)
    import numpy as np
    pts = m(np.array([0.1, 0.2]), np.array([0.1, 0.2]))
    assert gf(pts)[:, 0] == pytest.approx([0.01, 0.04])
    n = specialcf.normal(2)
    assert n(m(1, 0.5, VOL_or_BND=BND)) == pytest.approx((1, 0))


def test_deformation_and_replace(m):
    VH = VectorH1(m, order=2)
    d = GridFunction(VH)
    d.Set((0.1 * x, 0))
    m.SetDeformation(d)
    area = Integrate(1, m)
    m.UnsetDeformation()
    assert area == pytest.approx(1.1)
    assert Integrate(CF(1) * dx(deformation=d), m) == pytest.approx(1.1)
    H = H1(m, order=2)
    u = H.TrialFunction()
    gf = GridFunction(H)
    gf.Set(x)
    assert Integrate(grad(u).ReplaceFunction(gf)[0], m) == pytest.approx(1)
    assert Integrate(gf.CF(u.__diffop__()), m) == pytest.approx(0.5)


# ---------------------------------------------------------------- negative cases

@pytest.fixture(scope="module")
def spaces(m):
    H = H1(m, order=2)
    D = HDiv(m, order=1)
    L = L2(m, order=1, dgjumps=True)
    X = H1(m, order=1) * HDiv(m, order=1)
    gf = GridFunction(H)
    gf.Set(x)
    return H, D, L, X, gf


def test_hdiv_test_in_h1_linearform(spaces):
    H, D, L, X, gf = spaces
    t = D.TestFunction()
    f = LinearForm(H)
    raises(lambda: f.__iadd__(CF((1, 1)) * t * dx), "HDiv")


def test_other_family_in_forms(spaces):
    H, D, L, X, gf = spaces
    u, v = H.TnT()
    s, t = D.TnT()
    a = BilinearForm(D)
    raises(lambda: a.__iadd__(grad(u) * grad(v) * dx))
    b = BilinearForm(H)
    raises(lambda: b.__iadd__(Variation(s * s * dx)))
    c = BilinearForm(L)
    n = specialcf.normal(2)
    raises(lambda: c.__iadd__((s * n - s.Other() * n) * (t * n) * dx(skeleton=True)))
    d = BilinearForm(H)
    raises(lambda: d.__iadd__((u + CF((1, 1)) * s) * v * dx))
    e = BilinearForm(H)
    raises(lambda: e.__iadd__(InnerProduct(HDivDiv(H.mesh, order=1).TrialFunction(),
                                           HDivDiv(H.mesh, order=1).TestFunction()) * dx))
    F = BilinearForm(H)
    N = NumberSpace(H.mesh)
    nu, nv = N.TnT()
    raises(lambda: F.__iadd__(nu * nv * dx))


def test_compound_component_out_of_range(m):
    Y = H1(m, order=1) * H1(m, order=1) * L2(m) * L2(m)
    us, vs = Y.TnT()
    Z = H1(m, order=1) * H1(m, order=1)
    a = BilinearForm(Z)
    raises(lambda: a.__iadd__(us[3] * vs[3] * dx), "component")


def test_proxy_presence(spaces):
    H, D, L, X, gf = spaces
    u, v = H.TnT()
    a = BilinearForm(H)
    raises(lambda: a.__iadd__(u * dx), "test")
    raises(lambda: a.__iadd__(x * dx), "test")
    f = LinearForm(H)
    raises(lambda: f.__iadd__(u * v * dx), "rial")
    raises(lambda: f.__iadd__(x * dx), "test")
    raises(lambda: gf.Set(1 + u * u))
    raises(lambda: Integrate(u * v, H.mesh))
    raises(lambda: Integrate(u * v * dx, H.mesh))
    raises(lambda: u(H.mesh(0.2, 0.2)))
    b = BilinearForm(H)
    raises(lambda: b.__iadd__(Variation(u * v * dx)), "test")


def test_geometric_dimensions(m):
    raises(lambda: Integrate(specialcf.normal(3) * specialcf.normal(3), m, BND), "dimension")
    m3 = MakeStructured3DMesh(nx=2, ny=2, nz=2)
    raises(lambda: Integrate(specialcf.normal(2) * specialcf.normal(2), m3, BND), "dimension")
    raises(lambda: Integrate(specialcf.Weingarten(3), m, BND), "dimension")
    raises(lambda: Integrate(specialcf.tangential(3), m, BND), "dimension")
    raises(lambda: Integrate(specialcf.JacobianMatrix(2), m, BND), "dimension")


def test_domainwise_short_list(m2mat):
    raises(lambda: Integrate(CF([1]), m2mat), "regions")
    raises(lambda: Integrate(CF([7]), m2mat, BND), "regions")
    raises(lambda: Integrate(m2mat.BoundaryCF({".*": 1}), m2mat))


def test_complex_in_real_form(spaces):
    H, D, L, X, gf = spaces
    u, v = H.TnT()
    a = BilinearForm(H)
    raises(lambda: a.__iadd__(1j * u * v * dx), "complex")
    f = LinearForm(H)
    raises(lambda: f.__iadd__(1j * v * dx), "complex")
    raises(lambda: gf.Set(1j * x), "complex")


def test_facet_proxy_in_dx(m):
    F = FacetFESpace(m, order=1)
    uf, vf = F.TnT()
    a = BilinearForm(F)
    raises(lambda: a.__iadd__(uf * vf * dx))
    a += uf * vf * dx(element_boundary=True)
    a.Assemble()


def test_other_without_skeleton(spaces, m):
    H, D, L, X, gf = spaces
    u, v = L.TnT()
    a = BilinearForm(L)
    raises(lambda: a.__iadd__(SymbolicBFI(u.Other() * v)), "skeleton")
    raises(lambda: Integrate(x.Other(), m), "skeleton")


def test_dg_term_without_dgjumps(m):
    L = L2(m, order=1)
    u, v = L.TnT()
    a = BilinearForm(L)
    a += (u - u.Other()) * (v - v.Other()) * dx(skeleton=True)
    raises(lambda: a.Assemble(), "dgjumps")


def test_linearization_on_skeleton(spaces):
    H, D, L, X, gf = spaces
    u, v = L.TnT()
    a = BilinearForm(L)
    igl = (u - u.Other()) * v * dx(skeleton=True)
    igl.linearization = u * v * dx(skeleton=True)
    raises(lambda: a.__iadd__(igl), "linearization")


def test_deformation_checks(m):
    D = HDiv(m, order=1)
    gD = GridFunction(D)
    raises(lambda: m.SetDeformation(gD))
    m2 = Mesh(unit_square.GenerateMesh(maxh=0.4))
    g2 = GridFunction(VectorH1(m2, order=1))
    raises(lambda: m.SetDeformation(g2), "another mesh")
    gs = GridFunction(H1(m, order=1))
    raises(lambda: Integrate(CF(1) * dx(deformation=gs), m), "dim")


def test_replace_and_cf_diffop(spaces):
    H, D, L, X, gf = spaces
    s = D.TrialFunction()
    xs = X.TrialFunction()[1]
    raises(lambda: s.ReplaceFunction(gf))
    raises(lambda: div(s).ReplaceFunction(gf))
    raises(lambda: xs.ReplaceFunction(gf))
    raises(lambda: gf.CF(s.__diffop__()))


def test_interpolate_mixed_spaces(spaces):
    H, D, L, X, gf = spaces
    u = H.TrialFunction()
    s = D.TrialFunction()
    raises(lambda: Interpolate(u * CF((1, 1)) + s, VectorL2(H.mesh, order=0)), "different spaces")
