"""
Trial- and test-functions carry one differential operator per co-dimension of
the element they are evaluated on. In a form, Trace() is needed where the
restriction to a lower dimensional element is a projection (tangential trace of
HCurl, normal trace of HDiv, surface gradient of H1) and optional where it is
just the value (H1) or where the space lives on the boundary anyway (SurfaceL2,
HDivSurface, ...). GridFunctions and BEM potentials give what the element sees.
"""
import pickle
import pytest
from netgen.occ import OCCGeometry, Box
from ngsolve import *


@pytest.fixture(scope="module")
def mesh():
    return Mesh(OCCGeometry(Box((0, 0, 0), (1, 1, 1))).GenerateMesh(maxh=0.5))


def mat(fes, integrand):
    a = BilinearForm(fes, check_unused=False)
    a += integrand
    return a.Assemble().mat


def quad(m, gf):
    tmp = gf.vec.CreateVector()
    tmp.data = m*gf.vec
    return InnerProduct(tmp, gf.vec)


def same(m1, m2):
    diff = m1.AsVector().CreateVector()
    diff.data = m1.AsVector() - m2.AsVector()
    return Norm(diff) < 1e-12*Norm(m1.AsVector())


def raises(func, match="Trace"):
    with pytest.raises(Exception, match=match):
        func()


# unit cube, constant field (1,2,3)
field = CF((1, 2, 3))
tangential_on_faces = 2*(13+10+5)      # |v - (v.n)n|^2 over the 6 faces
normal_on_faces = 2*(1+4+9)            # (v.n)^2 over the 6 faces
tangential_on_edges = 4*(1+4+9)        # (v.t)^2 over the 12 edges


# ------------------------------------------------------------ H1: trace is the value

def test_h1_trace_is_optional(mesh):
    fes = H1(mesh, order=2)
    u, v = fes.TnT()
    gf = GridFunction(fes)
    gf.Set(1)
    m = mat(fes, u*v*ds)
    assert same(m, mat(fes, u.Trace()*v.Trace()*ds))
    assert quad(m, gf) == pytest.approx(6)
    medge = mat(fes, u*v*dx(definedon=mesh.BBoundaries(".*")))
    assert same(medge, mat(fes, u.Trace().Trace()*v.Trace().Trace()*dx(definedon=mesh.BBoundaries(".*"))))
    assert quad(medge, gf) == pytest.approx(12)


def test_h1_surface_gradient_needs_trace(mesh):
    fes = H1(mesh, order=2)
    u, v = fes.TnT()
    raises(lambda: mat(fes, grad(u)*grad(v)*ds))
    gf = GridFunction(fes)
    gf.Set(x+2*y+3*z)
    m = mat(fes, grad(u).Trace()*grad(v).Trace()*ds)
    assert quad(m, gf) == pytest.approx(tangential_on_faces)
    assert same(m, mat(fes, u.Trace().Deriv()*v.Trace().Deriv()*ds))


def test_h1_trace_chain(mesh):
    fes = H1(mesh, order=2)
    u = fes.TrialFunction()
    assert u.Trace().Trace().Trace() is not None
    assert "hesseboundary" in str(u.Operator("hesse").Trace())
    assert u.Deriv() is u.Deriv()


@pytest.mark.parametrize("make", [lambda m: H1(m, order=1),
                                  lambda m: H1(m, order=1, dim=2),
                                  lambda m: H1(m, order=1, dim=2, complex=True),
                                  lambda m: MatrixValued(H1(m, order=1), dim=3)])
def test_gridfunction_trace_keeps_dims(mesh, make):
    # the third trace lives on vertices (BBBND)
    gf = GridFunction(make(mesh))
    tr = gf
    for level in range(3):
        tr = tr.Trace()
        assert tuple(tr.dims) == tuple(gf.dims), level
    if gf.dim > 1:
        tr[gf.dim-1]


# ------------------------------------------------------------ HCurl, HDiv: trace is a projection

def test_hcurl_needs_trace(mesh):
    fes = HCurl(mesh, order=1)
    e, f = fes.TnT()
    gf = GridFunction(fes)
    gf.Set(field)
    raises(lambda: mat(fes, e*f*ds))
    assert quad(mat(fes, e.Trace()*f.Trace()*ds), gf) == pytest.approx(tangential_on_faces)

    edges = dx(definedon=mesh.BBoundaries(".*"))
    raises(lambda: mat(fes, e*f*edges))
    raises(lambda: mat(fes, e.Trace()*f.Trace()*edges))
    medge = mat(fes, e.Trace().Trace()*f.Trace().Trace()*edges)
    assert quad(medge, gf) == pytest.approx(tangential_on_edges)
    raises(lambda: e.Trace().Trace().Trace(), match="trace")

    raises(lambda: mat(fes, e.Trace()*f.Trace()*dx), match="lives on BND")
    raises(lambda: mat(fes, curl(e)*curl(f)*ds))
    assert same(mat(fes, curl(e).Trace()*curl(f).Trace()*ds),
                mat(fes, e.Trace().Deriv()*f.Trace().Deriv()*ds))


def test_hdiv_needs_trace(mesh):
    fes = HDiv(mesh, order=1)
    s, t = fes.TnT()
    gf = GridFunction(fes)
    gf.Set(field)
    raises(lambda: mat(fes, s*t*ds))
    assert quad(mat(fes, s.Trace()*t.Trace()*ds), gf) == pytest.approx(normal_on_faces)


def test_compound_components(mesh):
    fes = H1(mesh, order=2)*HCurl(mesh, order=1)
    (u, e), (v, f) = fes.TnT()
    gf = GridFunction(fes)
    gf.components[0].Set(1)
    gf.components[1].Set(field)
    raises(lambda: mat(fes, e*f*ds))
    m = mat(fes, (u*v + e.Trace()*f.Trace())*ds)
    assert quad(m, gf) == pytest.approx(6 + tangential_on_faces)


# ------------------------------------------------------------ GridFunctions give what the element sees

def test_gridfunction_on_boundary_is_the_trace(mesh):
    ge = GridFunction(HCurl(mesh, order=1))
    ge.Set(field)
    assert Integrate(ge*ge*ds, mesh) == pytest.approx(tangential_on_faces)
    assert Integrate(ge.Trace()*ge.Trace()*ds, mesh) == pytest.approx(tangential_on_faces)
    assert Integrate(ge*ge*dx, mesh) == pytest.approx(14)
    raises(lambda: Integrate(ge.Trace()*ge.Trace()*dx, mesh))

    gs = GridFunction(HDiv(mesh, order=1))
    gs.Set(field)
    assert Integrate(gs*gs*ds, mesh) == pytest.approx(normal_on_faces)


def test_gridfunction_trace_goes_one_level_down(mesh):
    edges = dx(definedon=mesh.BBoundaries(".*"))
    ge = GridFunction(HCurl(mesh, order=1))
    ge.Set(field)
    tr = ge.Trace()
    assert Integrate(tr*tr*ds, mesh) == pytest.approx(tangential_on_faces)
    assert Integrate(tr*tr*edges, mesh) == pytest.approx(tangential_on_edges)
    trtr = tr.Trace()
    assert Integrate(trtr*trtr*edges, mesh) == pytest.approx(tangential_on_edges)
    with pytest.raises(Exception):
        Integrate(trtr*trtr*ds, mesh)
    raises(lambda: trtr.Trace(), match="trace")
    raises(lambda: GridFunction(L2(mesh, order=1)).Trace(), match="trace")

    gh = GridFunction(H1(mesh, order=2))
    gh.Set(x+2*y+3*z)
    surfgrad = grad(gh).Trace()
    assert Integrate(surfgrad*surfgrad*ds, mesh) == pytest.approx(tangential_on_faces)
    hb = gh.Operator("hesseboundary")
    raises(lambda: hb.Trace(), match="trace")


def test_replacefunction_keeps_the_trace(mesh):
    fes = HCurl(mesh, order=1)
    e = fes.TrialFunction()
    ge = GridFunction(fes)
    ge.Set(field)
    full = e.ReplaceFunction(ge)
    assert Integrate(full*full*dx, mesh) == pytest.approx(14)
    assert Integrate(full*full*ds, mesh) == pytest.approx(tangential_on_faces)
    tr = e.Trace().ReplaceFunction(ge)
    assert Integrate(tr*tr*ds, mesh) == pytest.approx(tangential_on_faces)
    with pytest.raises(Exception):
        Integrate(tr*tr*dx, mesh)


def test_gridfunction_diff_through_trace(mesh):
    fes = H1(mesh, order=2)
    v = fes.TestFunction()
    gf = GridFunction(fes)
    gf.Set(x+2*y+3*z)

    def vec(integrand):
        f = LinearForm(fes)
        f += integrand
        return f.Assemble().vec

    ref = vec(gf*v*ds)
    d = vec((0.5*gf.Trace()*gf.Trace()).Diff(gf, v)*ds)
    d.data -= ref
    assert Norm(d) < 1e-12*Norm(ref)


def test_gridfunction_operator_pickle(mesh):
    gf = GridFunction(H1(mesh, order=2))
    gf.Set(x*x+y*z)
    # Hessian: xx = 2, yz = zy = 1; its tangential part on the 6 faces
    for op, d, val in [("hesse", dx, 6), ("hesseboundary", ds, 2*(2+4+4))]:
        cf = gf.Operator(op)
        assert Integrate(InnerProduct(cf, cf)*d, mesh) == pytest.approx(val)
        mesh2, cf2 = pickle.loads(pickle.dumps((mesh, cf)))
        assert Integrate(InnerProduct(cf2, cf2)*d, mesh2) == pytest.approx(val)


# ------------------------------------------------------------ surface spaces live on the boundary

def test_surface_l2(mesh):
    fes = SurfaceL2(mesh, order=1)
    a, b = fes.TnT()
    gf = GridFunction(fes)
    gf.Set(1, definedon=mesh.Boundaries(".*"))
    m = mat(fes, a*b*ds)
    assert same(m, mat(fes, a.Trace()*b.Trace()*ds))
    assert quad(m, gf) == pytest.approx(6)
    stiffness = mat(fes, grad(a)*grad(b)*ds)
    assert Norm(stiffness.AsVector()) > 0
    assert same(stiffness, mat(fes, grad(a).Trace()*grad(b).Trace()*ds))
    assert same(stiffness, mat(fes, a.Trace().Deriv()*b.Trace().Deriv()*ds))
    raises(lambda: mat(fes, a*b*dx), match="lives on BND")
    assert Integrate(gf*ds, mesh) == pytest.approx(6)
    assert Integrate(gf.Trace()*ds, mesh) == pytest.approx(6)
    # the GridFunction is zero inside the volume elements
    assert Integrate(gf*dx, mesh) == 0
    with pytest.raises(Exception, match="definedon"):
        gf.Set(1)


def test_hdiv_surface(mesh):
    fes = HDivSurface(mesh, order=1)
    s, t = fes.TnT()
    m = mat(fes, s*t*ds)
    assert Norm(m.AsVector()) > 0
    assert same(m, mat(fes, s.Trace()*t.Trace()*ds))
    d = mat(fes, div(s)*div(t)*ds)
    assert Norm(d.AsVector()) > 0
    assert same(d, mat(fes, div(s.Trace())*div(t.Trace())*ds))
    raises(lambda: mat(fes, s*t*dx), match="lives on BND")


@pytest.mark.parametrize("space", [VectorSurfaceL2, HDivDivSurface, FacetSurface, NormalFacetSurface])
def test_surface_spaces_have_no_volume_form(mesh, space):
    fes = space(mesh, order=1)
    a, b = fes.TnT()
    raises(lambda: mat(fes, InnerProduct(a, b)*dx), match="lives on BND")
    a.Trace()   # still allowed


def test_defined_on_region(mesh):
    from ngsolve.comp import DefinedOnRegion

    def count(cf):
        return [sum(DefinedOnRegion(cf, reg).Mask())
                for reg in (mesh.Materials(".*"), mesh.Boundaries(".*"), mesh.BBoundaries(".*"))]

    nmat, nbnd, nbbnd = len(mesh.GetMaterials()), len(mesh.GetBoundaries()), len(mesh.GetBBoundaries())
    gh = GridFunction(H1(mesh, order=2))
    gs = GridFunction(SurfaceL2(mesh, order=1))
    assert count(gh) == [nmat, nbnd, nbbnd]
    assert count(gh.Trace()) == [0, nbnd, nbbnd]
    assert count(GridFunction(L2(mesh, order=1))) == [nmat, 0, 0]
    assert count(gs) == [0, nbnd, 0]
    assert count(grad(gs)) == [0, nbnd, 0]
    assert count(gs*gh + x) == [0, nbnd, 0]
    assert count(x) == [nmat, nbnd, nbbnd]


def test_webgui_skips_the_volume_of_surface_functions(mesh):
    webgui = pytest.importorskip("ngsolve.webgui")
    regions = {VOL: VOL, BND: BND, BBND: BBND}
    gh = GridFunction(H1(mesh, order=2))
    gh.Set(x)
    gs = GridFunction(SurfaceL2(mesh, order=1))
    gs.Set(x, definedon=mesh.Boundaries(".*"))
    d = webgui.BuildRenderData(mesh, gh, regions=regions)
    assert d["draw_vol"] and "points3d" in d
    for cf in [gs, gh.Trace(), gs*gh]:
        d = webgui.BuildRenderData(mesh, cf, regions=regions)
        assert not d["draw_vol"] and "points3d" not in d
        assert d["funcmax"] == pytest.approx(1, abs=1e-8)


def test_webgui_other_mesh():
    # region numbers of the drawing mesh say nothing about the space's mesh
    from netgen.occ import Glue
    from ngsolve.comp import DefinedOnRegion
    webgui = pytest.importorskip("ngsolve.webgui")
    left, right = Box((0, 0, 0), (0.5, 1, 1)), Box((0.5, 0, 0), (1, 1, 1))
    left.mat("left")
    right.mat("right")
    source = Mesh(OCCGeometry(Glue([left, right])).GenerateMesh(maxh=0.6))
    target = Mesh(OCCGeometry(Box((0, 0, 0), (1, 1, 1))).GenerateMesh(maxh=0.6))
    gf = GridFunction(H1(source, order=1, definedon=source.Materials("right")))
    gf.Set(1)
    assert gf(target(0.75, 0.5, 0.5)) == pytest.approx(1)
    assert list(DefinedOnRegion(gf, target.Materials(".*")).Mask()) == [True]
    d = webgui.BuildRenderData(target, gf, regions={VOL: VOL, BND: BND, BBND: BBND})
    assert d["draw_vol"] and "points3d" in d


def test_defined_on_region_elements_without_region():
    import netgen.meshing as nm
    from ngsolve.comp import DefinedOnRegion
    raw = nm.Mesh(3)
    pts = [raw.Add(nm.MeshPoint(nm.Pnt(*p))) for p in [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]]
    raw.Add(nm.Element3D(1, pts))
    face = raw.Add(nm.FaceDescriptor(surfnr=1, domin=1, bc=1))
    for t in [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]:
        raw.Add(nm.Element2D(face, [pts[i] for i in t]))
    # the first segment has no region (index -1)
    raw.Add(nm.Element1D([pts[0], pts[1]], index=0))
    raw.Add(nm.Element1D([pts[2], pts[3]], index=1))
    m = Mesh(raw)
    assert [el.index for el in m.Elements(BBND)] == [-1, 0]
    assert list(DefinedOnRegion(CF(1), m.BBoundaries(".*")).Mask()) == [True]


# ------------------------------------------------------------ BEM potentials take the trace themselves

@pytest.mark.parametrize("use_fmm", [False, True])
@pytest.mark.parametrize("untraced", ["trial", "test", "both"])
def test_bem_falls_back_to_the_trace(mesh, use_fmm, untraced):
    bem = pytest.importorskip("ngsolve.bem")
    fes = HCurl(mesh, order=1, complex=True)
    e, f = fes.TnT()
    gf = GridFunction(fes)
    gf.Set(field)

    def apply(trial, test):
        y = gf.vec.CreateVector()
        y.data = (bem.HelmholtzSL(trial*ds, 1.0, use_fmm=use_fmm)*test*ds).mat * gf.vec
        return y

    ref = apply(e.Trace(), f.Trace())
    y = apply(e if untraced != "test" else e.Trace(), f if untraced != "trial" else f.Trace())
    assert Norm(ref) > 0
    y.data -= ref
    assert Norm(y) < 1e-10*Norm(ref)


# ------------------------------------------------------------ ConvertOperator takes the operator of the elements

def test_convertoperator_with_traced_proxy(mesh):
    fes = HCurl(mesh, order=1)
    e = fes.TrialFunction()
    gf = GridFunction(fes)
    gf.Set(field)

    def apply(**kwargs):
        y = gf.vec.CreateVector()
        y.data = ConvertOperator(fes, fes, **kwargs) * gf.vec
        return y

    ref = apply(vb=BND)
    assert Norm(ref) > 0
    for proxy in [e, e.Trace()]:
        y = apply(trial_proxy=proxy, vb=BND)
        y.data -= ref
        assert Norm(y) < 1e-12*Norm(ref)
    with pytest.raises(Exception):
        apply(trial_proxy=e.Trace(), vb=VOL)


# ------------------------------------------------------------ the trace knows where it comes from

def test_trace_is_cached(mesh):
    u = H1(mesh, order=2).TrialFunction()
    assert u.Trace() is u.Trace()
    assert u.Trace().Trace() is u.Trace().Trace()
    assert u.Trace().primary is u
    assert grad(u).Trace().primary is grad(u)
    assert u.dt.Trace() is u.Trace().dt
    assert "order-dt=1" in str(u.dt.Trace())


def test_diff_through_trace(mesh):
    fes = H1(mesh, order=2)
    u, v = fes.TnT()
    mass = mat(fes, u*v*ds)
    assert same(mat(fes, (0.5*u.Trace()*u.Trace()).Diff(u, v)*ds), mass)
    assert same(mat(fes, (0.5*u.Trace()*u.Trace()).Diff(u.Trace(), v.Trace())*ds), mass)
    stiffness = mat(fes, grad(u).Trace()*grad(v).Trace()*ds)
    assert same(mat(fes, (0.5*grad(u).Trace()*grad(u).Trace()).Diff(u, v)*ds), stiffness)
    assert same(mat(fes, (0.5*u.Trace().Deriv()*u.Trace().Deriv()).Diff(u, v)*ds), stiffness)

    fes = H1(mesh, order=2)*HCurl(mesh, order=1)
    (u, e), (v, f) = fes.TnT()
    energy = 0.5*u*u + 0.5*e.Trace()*e.Trace()
    assert same(mat(fes, (energy.Diff(u, v) + energy.Diff(e, f))*ds),
                mat(fes, (u*v + e.Trace()*f.Trace())*ds))
