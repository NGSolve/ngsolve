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


# ---------------------------------------------------------------- stage 2: more entry points

def test_numpy_point_arrays(m):
    u = H1(m, order=1).TrialFunction()
    pts = m.MapToAllElements(IntegrationRule(TRIG, 2), VOL)
    raises(lambda: u(pts), "trial")
    raises(lambda: specialcf.normal(3)(pts), "dimension")
    raises(lambda: specialcf.Weingarten(2)(pts), "facet")
    assert (x * y)(pts).shape == (len(pts), 1)


def test_vtk_output(m, tmp_path):
    u = H1(m, order=1).TrialFunction()
    for cf in [u, specialcf.normal(3), CF([1, 2, 3])[0] * 0 + m.BoundaryCF({".*": 1})]:
        vtk = VTKOutput(m, coefs=[cf], names=["f"], filename=str(tmp_path / "out"))
        raises(lambda: vtk.Do())
    vtk = VTKOutput(m, coefs=[x * y], names=["f"], filename=str(tmp_path / "ok"))
    vtk.Do()
    vtk.Do(BND)


def test_mass_operators(m, m2mat):
    L = L2(m, order=1)
    u = L.TrialFunction()
    vec = GridFunction(L).vec
    raises(lambda: L.Mass(rho=u))
    raises(lambda: L.SolveM(vec, rho=specialcf.normal(3)), "dimension")
    raises(lambda: L.ApplyM(vec, rho=1j * x), "complex")
    raises(lambda: L2(m2mat, order=1).Mass(rho=CF([1])), "regions")
    gh = GridFunction(H1(m, order=1))
    M = gh.space.Mass(1)
    raises(lambda: gh.vec.CreateVector().data.__class__ and M.Mult(gh.vec, gh.vec.CreateVector()), "L2")
    L.SolveM(vec, rho=1 + x)


def test_classic_integrator_coefficients(m, m2mat):
    H = H1(m, order=1)
    u = H.TrialFunction()
    a = BilinearForm(H)
    raises(lambda: a.__iadd__(BFI("mass", coef=u)))
    a2 = BilinearForm(H1(m2mat, order=1))
    raises(lambda: a2.__iadd__(BFI("mass", coef=CF([1]))), "regions")
    f = LinearForm(H)
    raises(lambda: f.__iadd__(LFI("source", coef=specialcf.normal(3)[0])), "dimension")
    a += BFI("mass", coef=1 + x)
    a.Assemble()


def test_placeholder_set(m):
    H = H1(m, order=1)
    u, v = H.TnT()
    p = PlaceholderCF(x)
    assert Integrate(p, m) == pytest.approx(0.5)
    assert p(m(0.3, 0.4)) == pytest.approx(0.3)
    raises(lambda: p.Set(u), "trial")
    raises(lambda: p.Set(specialcf.normal(3)[0]), "dimension")
    raises(lambda: p.Set(1j * x), "real")
    f = LinearForm(H)
    q = PlaceholderCF(x)
    f += q * v * dx
    raises(lambda: q.Set(u), "trial")
    q.Set(y)
    f.Assemble()
    p.Set(y)
    assert p(m(0.3, 0.4)) == pytest.approx(0.4)


def test_newton_cf_expression(m):
    from ngsolve.fem import NewtonCF
    L = L2(m, order=0)
    u = L.TrialFunction()
    ncf = NewtonCF(u * u - 2, 1)
    assert Integrate(ncf, m) == pytest.approx(2 ** 0.5)
    bad = NewtonCF(u * u - 2 + x * specialcf.normal(3)[0], 1)
    raises(lambda: Integrate(bad, m), "dimension")
    v = H1(m, order=1).TrialFunction()
    raises(lambda: Integrate(NewtonCF(u * u - 2, v), m), "trial")


def test_interpolate_cf_context(m):
    V = H1(m, order=2)
    raises(lambda: Integrate(Interpolate(x, V), m, BND), "Interpolate")
    raises(lambda: Integrate(Interpolate(specialcf.normal(3)[0], V), m), "dimension")
    assert Integrate(Interpolate(x * y, V), m) == pytest.approx(0.25)


def test_pml_custom(m):
    u = H1(m, order=1).TrialFunction()
    raises(lambda: m.SetPML(pml.Custom(trafo=CF((x + 1j * u, y)), jac=CF((1, 0, 0, 1), dims=(2, 2))), 1))


def test_convert_operator(m):
    F = FacetFESpace(m, order=1)
    L = L2(m, order=1)
    raises(lambda: ConvertOperator(F, L), "inside elements")
    H = H1(m, order=1)
    raises(lambda: ConvertOperator(H, L, trial_cf=H.TrialFunction() * specialcf.normal(3)[0]), "dimension")
    ConvertOperator(H, L)


def test_point_evaluation_functional(m):
    H = H1(m, order=1)
    u, v = H.TnT()
    f = LinearForm(H)
    raises(lambda: f.__iadd__((u * v)(0.5, 0.5)), "rial")
    raises(lambda: f.__iadd__(CF((v, v))(0.5, 0.5)), "scalar")
    raises(lambda: f.__iadd__(x(0.5, 0.5)), "test")
    f2 = LinearForm(H)
    f2 += v(5, 5)
    raises(lambda: f2.Assemble(), "not in the mesh")


def test_facet_spaces_inside_elements(m):
    for fes in [FacetFESpace(m, order=1), TangentialFacetFESpace(m, order=1), NormalFacetFESpace(m, order=1)]:
        gf = GridFunction(fes)
        gf.Set(x if fes.type == "facet" else CF((x, y)), dual=True)
        raises(lambda: Integrate(gf, m))
        raises(lambda: gf(m(0.3, 0.4)))
        u, v = fes.TnT()
        a = BilinearForm(fes)
        raises(lambda: a.__iadd__(u * v * dx), "inside elements")
        Integrate(InnerProduct(gf, gf) * dx(element_boundary=True), m)
    F = FacetFESpace(m, order=1)
    gF = GridFunction(F)
    raises(lambda: gF.Set(x), "dual")
    raises(lambda: gF(0.3, 0.4), "inside")


def test_facet_surface_on_boundary():
    m3 = MakeStructured3DMesh(nx=2, ny=2, nz=2)
    F = FacetSurface(m3, order=1)
    u, v = F.TnT()
    a = BilinearForm(F)
    raises(lambda: a.__iadd__(u * v * ds), "inside elements")
    a += u * v * ds(element_boundary=True)
    a.Assemble()


def test_weingarten_and_tangent_need_facets(m):
    raises(lambda: Integrate(specialcf.Weingarten(2)[0, 0], m), "facet")
    Integrate(specialcf.Weingarten(2)[0, 0] * dx(element_boundary=True), m)
    raises(lambda: Integrate(specialcf.tangential(2, consistent=True)[0], m), "facet")


def test_gridfunction_after_refinement():
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.4))
    H = H1(mesh, order=1)
    gf = GridFunction(H, autoupdate=False)
    gf.Set(x)
    u, v = H.TnT()
    a = BilinearForm(H)
    a += gf * u * v * dx
    f = LinearForm(H)
    f += gf * v * dx
    cf = gf * x
    cf(mesh(0.3, 0.3))
    mesh.Refine()
    H.Update()
    raises(lambda: a.Assemble(), "up to date")
    raises(lambda: f.Assemble(), "up to date")
    raises(lambda: cf(mesh(0.3, 0.3)), "up to date")
    raises(lambda: gf(0.3, 0.3), "up to date")
    gf.Update()
    a.Assemble()
    f.Assemble()
    cf(mesh(0.3, 0.3))


def test_exception_in_assembly_stops(m):
    from ngsolve.fem import NewtonCF
    L = L2(m, order=0)
    w = L.TrialFunction()
    ncf = NewtonCF(w * w + 1, 1, maxiter=2)
    H = H1(m, order=1)
    u, v = H.TnT()
    a = BilinearForm(H)
    a += ncf * u * v * dx
    raises(lambda: a.Assemble())


def test_proxy_names_demangled(spaces):
    H, D, L, X, gf = spaces
    t = D.TestFunction()
    f = LinearForm(H)
    with pytest.raises(Exception) as e:
        f.__iadd__(CF((1, 1)) * t * dx)
    assert "FN5ngfem" not in str(e.value)


def test_dof_consistency(m):
    # the assembly loops rely on GetFE and GetDofNrs agreeing, checked once per space update
    m3 = MakeStructured3DMesh(nx=2, ny=2, nz=2)
    for mesh in [m, m3]:
        for fes in [H1(mesh, order=3), L2(mesh, order=2), HDiv(mesh, order=2), HCurl(mesh, order=2),
                    FacetFESpace(mesh, order=2), HDivDiv(mesh, order=1), VectorH1(mesh, order=2),
                    H1(mesh, order=2) * HDiv(mesh, order=1), Discontinuous(H1(mesh, order=2))]:
            fes._CheckDofConsistency(VOL)
            fes._CheckDofConsistency(BND)


# ------------------------------------------------- every CF node reports its inputs

def _cf_nodes():
    from ngsolve.fem import Einsum, CoordinateTrafo, PlaceholderCF
    from ngsolve.comp import BoundaryFromVolumeCF
    nodes = {
        "add": lambda u, M, s, m: u + u,
        "scale": lambda u, M, s, m: 2 * u,
        "scale_complex": lambda u, M, s, m: 2j * s,
        "mult_scal_vec": lambda u, M, s, m: s * u,
        "inner": lambda u, M, s, m: InnerProduct(u, u),
        "norm": lambda u, M, s, m: Norm(u),
        "norm_complex": lambda u, M, s, m: Norm(1j * u),
        "mat_vec": lambda u, M, s, m: M * u,
        "mat_mat": lambda u, M, s, m: M * M,
        "cross": lambda u, M, s, m: Cross(CF((s, s, s)), CF((1, 2, 3))),
        "component": lambda u, M, s, m: u[0],
        "subtensor": lambda u, M, s, m: M[0, :],
        "transpose": lambda u, M, s, m: M.trans,
        "trace": lambda u, M, s, m: Trace(M),
        "det": lambda u, M, s, m: Det(M),
        "inv": lambda u, M, s, m: Inv(M),
        "cof": lambda u, M, s, m: Cof(M),
        "sym": lambda u, M, s, m: Sym(M),
        "skew": lambda u, M, s, m: Skew(M),
        "identity_mult": lambda u, M, s, m: Id(2) * s,
        "unary": lambda u, M, s, m: sin(s),
        "binary": lambda u, M, s, m: atan2(s, 1),
        "power": lambda u, M, s, m: s ** 2,
        "ifpos": lambda u, M, s, m: IfPos(s, u, -u),
        "ifpos_condition": lambda u, M, s, m: IfPos(s, 1, 2),
        "vectorial": lambda u, M, s, m: CF((s, 2 * s)),
        "reshape": lambda u, M, s, m: u.Reshape((2, 1)),
        "tensor_transpose": lambda u, M, s, m: M.TensorTranspose((1, 0)),
        "extend_dimension": lambda u, M, s, m: u.ExtendDimension((3,)),
        "einsum": lambda u, M, s, m: Einsum("ij,j->i", M, u),
        "compile": lambda u, M, s, m: (s * s).Compile(),
        "eig": lambda u, M, s, m: M.Eig(),
        "boundary_from_volume": lambda u, M, s, m: BoundaryFromVolumeCF(u),
        "freeze": lambda u, M, s, m: s.Freeze(),
        "real": lambda u, M, s, m: s.real,
        "imag": lambda u, M, s, m: s.imag,
        "conj": lambda u, M, s, m: Conj(s),
        "parameter": lambda u, M, s, m: Parameter(2) * s,
        "domainwise": lambda u, M, s, m: m.MaterialCF({".*": s}),
        "regionwise": lambda u, M, s, m: m.RegionCF(VOL, {".*": s}),
        "coordinate_trafo": lambda u, M, s, m: s(CoordinateTrafo(CF((x, y)), m.Materials(".*"))),
        "placeholder": lambda u, M, s, m: PlaceholderCF(s),
        "bspline": lambda u, M, s, m: BSpline(2, [0, 0, 1, 2, 3], [1, 2, 3, 4])(s),
        "other": lambda u, M, s, m: u.Other(),
        "diff": lambda u, M, s, m: (s * s * x).Diff(x),
    }
    return nodes


@pytest.mark.parametrize("name", sorted(_cf_nodes().keys()))
def test_cf_node_reports_proxies(name):
    m = MakeStructured2DMesh(nx=2, ny=2)
    u = VectorH1(m, order=2).TrialFunction()
    node = _cf_nodes()[name](u, grad(u), u[0], m)
    msg = "can only be used in forms"
    raises(lambda: Integrate(node, m), msg)
    if name not in ("domainwise", "regionwise"):
        raises(lambda: Integrate(node, m, BND), msg)
    raises(lambda: node(m(0.3, 0.3)), msg)


def test_hidden_inputs_now_visible(m):
    # finding: BoundaryFromVolumeCF and Eig hid proxies from the gate
    from ngsolve.comp import BoundaryFromVolumeCF
    V = H1(m, order=2)
    u, v = V.TnT()
    a1 = BilinearForm(V)
    a1 += BoundaryFromVolumeCF(u) * v * dx
    a1.Assemble()
    a2 = BilinearForm(u * v * dx).Assemble()
    assert abs(nrm(a1) - nrm(a2)) < 1e-10 * nrm(a2)
    raises(lambda: BilinearForm(V).__iadd__(BoundaryFromVolumeCF(u) * v * ds), "can only be used in forms")
    gf = GridFunction(L2(m, order=2))
    gf.Set(x * x * y)
    bfv = BoundaryFromVolumeCF(grad(gf))[0]
    i1 = Integrate(bfv * ds, m)
    assert abs(i1 - Integrate(bfv.Compile() * ds, m)) < 1e-12
    M = CF((1 + x * x, x, x, 2), dims=(2, 2))
    gM = GridFunction(MatrixValued(L2(m, order=2)))
    gM.Set(M)
    eig = gM.Eig()
    Integrate(eig, m)
    p = Parameter(1)
    assert Integrate((p * eig[4]).Diff(p) - eig[4], m) == pytest.approx(0, abs=1e-12)
    raises(lambda: Integrate(CF((1j * x, 0, 0, 1), dims=(2, 2)).Eig(), m), "real")


def test_hidden_inputs_replace(m):
    from ngsolve.comp import BoundaryFromVolumeCF
    gf = GridFunction(H1(m, order=1))
    gf.Set(x)
    bfv = BoundaryFromVolumeCF(gf) * x
    r = bfv.Replace({x: y}) if hasattr(bfv, "Replace") else bfv
    Integrate(r * ds, m)


# ------------------------------------------------- operators of another mesh dimension

def test_proxy_of_other_dimension():
    m2 = MakeStructured2DMesh(nx=3, ny=3)
    fes = H1(m2, order=2)
    u, v = fes.TnT()
    u3 = H1(MakeStructured3DMesh(hexes=False, nx=1, ny=1, nz=1), order=2).TrialFunction()
    raises(lambda: BilinearForm(fes).__iadd__(u3 * v * dx), "3D space")
    raises(lambda: BilinearForm(fes).__iadd__(grad(u3)[0] * grad(v)[0] * dx), "3D space")


def test_proxy_of_other_mesh_same_dimension():
    ma = MakeStructured2DMesh(nx=3, ny=3)
    mb = MakeStructured2DMesh(nx=2, ny=2)
    fa = H1(ma, order=2)
    u, v = fa.TnT()
    ub = H1(mb, order=2).TrialFunction()
    a = BilinearForm(fa)
    a += grad(ub) * grad(v) * dx + ub * v * dx
    a.Assemble()
    b = BilinearForm(grad(u) * grad(v) * dx + u * v * dx).Assemble()
    assert abs(nrm(a) - nrm(b)) < 1e-10 * nrm(b)


def test_gridfunction_of_other_dimension():
    m2 = MakeStructured2DMesh(nx=3, ny=3)
    m3 = MakeStructured3DMesh(hexes=False, nx=2, ny=2, nz=2)
    g3 = GridFunction(H1(m3, order=1))
    g3.Set(x + y + z)
    raises(lambda: Integrate(g3, m2), "points of dimension")
    raises(lambda: g3(m2(0.3, 0.3)), "points of dimension")
    g2 = GridFunction(H1(m2, order=1))
    g2.Set(x + y)
    raises(lambda: Integrate(g2, m3), "points of dimension")
    # same dimension: point search
    mb = MakeStructured2DMesh(nx=5, ny=5)
    assert Integrate(g2, mb) == pytest.approx(1.0)
    # surface mesh (no volume elements) evaluated at points of a 3D mesh
    from ngsolve.meshes import MakeStructuredSurfaceMesh
    ms = MakeStructuredSurfaceMesh(nx=2, ny=2)
    gs = GridFunction(VectorL2(ms, order=1))
    mp = MakeStructured3DMesh(hexes=False, prism=True, nx=1, ny=1, nz=1)
    raises(lambda: gs(mp(0.0, 0.0, 0.0)), "another mesh")


# ------------------------------------------------- Other(), facet operators, element_vb

@pytest.mark.parametrize("space", ["h1", "l2", "discontinuous"])
def test_other_in_linearform_needs_skeleton(space):
    mesh = MakeStructured2DMesh(nx=2, ny=2)
    fes = {"h1": lambda: H1(mesh, order=1), "l2": lambda: L2(mesh, order=1),
           "discontinuous": lambda: Discontinuous(H1(mesh, order=1))}[space]()
    v = fes.TestFunction()
    for dX in [ds(element_vb=BND), dx(element_vb=BND), ds(element_boundary=True), dx(element_boundary=True)]:
        raises(lambda: LinearForm(y * v.Other() * dX), "Other")
    LinearForm(y * v.Other() * dx(skeleton=True)).Assemble()


def test_other_element_boundary_bilinear(m):
    fes = L2(m, order=1, dgjumps=True)
    u, v = fes.TnT()
    n = specialcf.normal(2)
    a = BilinearForm(fes)
    a += 0.5 * (u + u.Other()) * n[0] * v * dx(element_boundary=True)
    a.Assemble()
    raises(lambda: BilinearForm(fes).__iadd__(u.Other() * v * ds(element_boundary=True)), "Other")
    raises(lambda: Integrate(x.Other() * ds(element_boundary=True), m), "Other")
    Integrate(x.Other() * dx(element_boundary=True), m)


def test_facet_space_on_edges():
    m3 = MakeStructured3DMesh(hexes=True, nx=1, ny=1, nz=1)
    v = FacetFESpace(m3, order=1).TestFunction()
    raises(lambda: LinearForm(z * v * dx(element_vb=BBND)), "facets")
    LinearForm(z * v * dx(element_vb=BND)).Assemble()
    gf = GridFunction(FacetFESpace(m3, order=1))
    raises(lambda: Integrate(gf * dx(element_vb=BBND), m3), "inside elements")


def test_element_vb_sub_entities(m):
    # vertices of 2D elements are fine (dual bases), there is nothing below
    V = H1(m, order=1)
    u, v = V.TnT()
    BilinearForm(u * v * dx(element_vb=BBND)).Assemble()
    assert Integrate(CF(1) * dx(element_vb=BBND), m) == pytest.approx(3 * m.ne)
    raises(lambda: BilinearForm(V).__iadd__(u * v * dx(element_vb=BBBND)), "sub-entities")
    raises(lambda: BilinearForm(V).__iadd__(u * v * ds(element_vb=BBND)), "sub-entities")
    raises(lambda: Integrate(CF(1) * dx(element_vb=BBBND), m), "sub-entities")


# ------------------------------------------------- skeleton of surface meshes

def test_surface_mesh_boundary_skeleton():
    from ngsolve.meshes import MakeStructuredSurfaceMesh
    ms = MakeStructuredSurfaceMesh(nx=2, ny=2)
    u, v = H1(ms, order=1).TnT()
    raises(lambda: BilinearForm(u * v * ds(skeleton=True)), "volume elements")
    raises(lambda: LinearForm(v * ds(skeleton=True)), "volume elements")
    BilinearForm(u * v * ds).Assemble()


def test_boundary_skeleton_dg(m):
    fes = L2(m, order=1)
    u, v = fes.TnT()
    n = specialcf.normal(2)
    a = BilinearForm(fes)
    a += 10 * u * v * ds(skeleton=True) - n * grad(u) * v * ds(skeleton=True)
    a.Assemble()
    f = LinearForm(v * ds(skeleton=True)).Assemble()
    assert f.vec.Norm() > 0
    assert Integrate(CF(1) * ds(skeleton=True), m) == pytest.approx(4)


# ------------------------------------------------- geometric CFs

def test_num_els_on_facet():
    geo = g2d.SplineGeometry()
    geo.AddRectangle((0, 0), (1, 1), bcs=("bottom", "right", "top", "left"))
    mq = Mesh(geo.GenerateMesh(maxh=0.4, quad_dominated=True))
    ne = specialcf.num_els_on_facet()
    assert Integrate(ne * (1 + x), mq, BND) == pytest.approx(6)
    # interior facets count twice, boundary facets once: 2*interior + boundary length
    Integrate(ne * dx(element_boundary=True), mq)
    raises(lambda: Integrate(ne, mq), "facets")
    raises(lambda: Integrate(ne * dx(element_vb=BBND), mq), "facets")


def test_edge_face_tangential_vectors():
    m3 = MakeStructured3DMesh(hexes=False, nx=1, ny=1, nz=1)
    t = specialcf.EdgeFaceTangentialVectors(3)
    Integrate(t[0, 0] * dx(element_vb=BBND), m3)
    raises(lambda: Integrate(t[0, 0], m3), "edges")
    mh = MakeStructured3DMesh(hexes=True, nx=1, ny=1, nz=1)
    raises(lambda: Integrate(t[0, 0] * dx(element_vb=BBND), mh), "tets")


# ------------------------------------------------- definedonelements sizes

def test_definedonelements_sizes(m):
    from pyngcore import BitArray
    V = H1(m, order=1)
    u, v = V.TnT()
    short = BitArray(3)
    short.Set()
    raises(lambda: Integrate(CF(1) * dx(definedonelements=short), m), "definedonelements")
    raises(lambda: Integrate(x * dx(definedonelements=BitArray(0)), m), "definedonelements")
    raises(lambda: BilinearForm(V).__iadd__(u * v * dx(definedonelements=short)), "definedonelements")
    raises(lambda: LinearForm(V).__iadd__(v * ds(definedonelements=short)), "definedonelements")
    gf = GridFunction(V)
    raises(lambda: gf.Set(x, definedonelements=short))
    L = L2(m, order=1, dgjumps=True)
    p, q = L.TnT()
    ba = BitArray(m.ne)
    ba.Set()
    raises(lambda: BilinearForm(L).__iadd__(p * q * dx(skeleton=True, definedonelements=ba)), "facets")
    bf = BitArray(m.nfacet)
    bf.Set()
    BilinearForm(L).__iadd__(p * q * dx(skeleton=True, definedonelements=bf))
    BilinearForm(L).__iadd__(p * q * dx(element_boundary=True, definedonelements=ba))
    full = BitArray(m.ne)
    full.Set()
    assert Integrate(CF(1) * dx(definedonelements=full), m) == pytest.approx(1)


def test_other_boundary_values(m):
    fes = L2(m, order=1)
    u, v = fes.TnT()
    a = BilinearForm(fes, nonassemble=True)
    a += (u.Other(bnd=x + y) - u) * v * dx(skeleton=True) + (u.Other(bnd=x + y) - u) * v * ds(skeleton=True)
    gf = GridFunction(fes)
    gf.Set(x)
    r = gf.vec.CreateVector()
    a.Apply(gf.vec, r)
    raises(lambda: BilinearForm(fes).__iadd__(u.Other(bnd=v) * v * ds(skeleton=True)), "can only be used in forms")
    w = L2(m, order=0).TrialFunction()
    raises(lambda: BilinearForm(fes).__iadd__(u.Other(bnd=w) * v * ds(skeleton=True)), "can only be used in forms")
