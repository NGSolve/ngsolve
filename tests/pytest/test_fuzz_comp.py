"""
Fuzzing regressions in comp: illegal constructs raise clean exceptions,
1D variants of spaces work, results stay correct.
"""
import pickle
import pytest
import ngsolve
from ngsolve import *
from ngsolve.comp import Hidden, Reorder
from ngsolve.meshes import Make1DMesh, MakeStructured2DMesh, MakeStructured3DMesh, MakeStructuredSurfaceMesh
from netgen.geom2d import unit_square
from netgen.csg import unit_cube


def l2err(gf, f, mesh, vb=VOL):
    return sqrt(Integrate(InnerProduct(gf-f, gf-f), mesh, vb))


# ------------------------------------------------------------ spaces in 1D

def test_h1_1d_dim():
    m = Make1DMesh(4)
    fes = H1(m, order=2, dim=2)
    gf = GridFunction(fes)
    gf.Set((x, x*x))
    assert l2err(gf, CF((x, x*x)), m) < 1e-12
    assert l2err(Grad(gf), CF((1, 2*x), dims=(2,1)), m) < 1e-12


def test_vectorl2_1d():
    m = Make1DMesh(3)
    fes = VectorL2(m, order=2)
    gf = GridFunction(fes)
    gf.Set(CF((x*x,)))
    assert l2err(gf, CF((x*x,)), m) < 1e-12
    assert l2err(div(gf), 2*x, m) < 1e-12
    GridFunction(fes).Operator("dual")
    u, v = fes.TnT()
    assert tuple(u.dims) == (1,)
    VectorValued(Hidden(fes), dim=2).TnT()
    MatrixValued(fes, 2).TnT()
    with pytest.raises(Exception):
        VectorL2(m, order=1, piola=True)


def test_vectorh1_1d():
    m = Make1DMesh(2)
    fes = VectorH1(m, order=2)
    gf = GridFunction(fes)
    gf.Set(CF((x*x,)))
    assert l2err(gf, CF((x*x,)), m) < 1e-12
    assert l2err(Grad(gf), CF((2*x,), dims=(1,1)), m) < 1e-12
    d = gf.Diff(gf, gf)
    assert l2err(d, CF((x*x,)), m) < 1e-12


def test_interpolate_facet_1d_no_grad():
    m = Make1DMesh(3)
    v = H1(m, order=1).TestFunction()
    with pytest.raises(Exception):
        Grad(Interpolate(v, FacetFESpace(m)))


def test_facet_1d_highest_order_dc():
    m = Make1DMesh(3)
    assert FacetFESpace(m, order=1, highest_order_dc=True).ndof == FacetFESpace(m, order=1).ndof


# ------------------------------------------------------------ NodalFESpace

@pytest.mark.parametrize("order", [0, 3])
def test_nodal_illegal_order(order):
    m = MakeStructured2DMesh(quads=False, nx=2, ny=2, secondorder=True)
    with pytest.raises(Exception):
        NodalFESpace(m, order=order)


def nodal_meshes():
    return [("segm", lambda: Make1DMesh(4)),
            ("trig", lambda: Mesh(unit_square.GenerateMesh(maxh=0.4))),
            ("trig2", lambda: MakeStructured2DMesh(quads=False, nx=2, ny=2, secondorder=True)),
            ("quad", lambda: MakeStructured2DMesh(quads=True, nx=3, ny=3)),
            ("tet", lambda: Mesh(unit_cube.GenerateMesh(maxh=0.5))),
            ("hex", lambda: MakeStructured3DMesh(hexes=True, nx=2, ny=2, nz=2))]


@pytest.mark.parametrize("name,mk", nodal_meshes())
def test_nodal_order2_exact(name, mk):
    m = mk()
    f = {1: x*x, 2: x*x+y-x*y, 3: x*x+y*z-x*y+z*z}[m.dim]
    fes = NodalFESpace(m, order=2)
    assert fes.ndof == m.nv + m.nedge
    g = GridFunction(fes)
    g.Set(f)
    assert l2err(g, f, m) < 1e-10
    if m.dim > 1:
        assert l2err(g, f, m, BND) < 1e-10


# ------------------------------------------------------------ other spaces

def test_hct_order():
    m = MakeStructured2DMesh(quads=False, nx=2, ny=2)
    for order in [0, 4]:
        with pytest.raises(Exception):
            HCT_FESpace(m, order=order)
    g = GridFunction(HCT_FESpace(m, order=3))
    g.Set(x*x*y)
    assert l2err(g, x*x*y, m) < 1e-10


def test_divfree_reconstruction():
    mt = MakeStructured2DMesh(quads=False, nx=2, ny=2)
    gf = GridFunction(VectorH1(mt, order=1))
    gf.Set((x, y))
    assert l2err(gf.Operator("divfree_reconstruction"), CF((x, y)), mt) < 1e-12
    gf = GridFunction(VectorH1(mt, order=2))
    gf.Set((y, x))
    assert l2err(gf.Operator("divfree_reconstruction"), CF((y, x)), mt) < 1e-12
    mq = MakeStructured2DMesh(quads=True, nx=2, ny=2)
    with pytest.raises(Exception):
        GridFunction(VectorH1(mq, order=2)).Operator("divfree_reconstruction")


def test_hdivdivfacet_2d_no_trace():
    m = MakeStructured2DMesh(quads=False, nx=2, ny=2)
    u, v = HDivDivFacetSpace(m, order=1).TnT()
    with pytest.raises(Exception):
        BilinearForm(InnerProduct(u.Trace(), v.Trace())*ds).Assemble()


def test_getorder_setorder_face_without_surface_element():
    fes = HDiv(MakeStructured2DMesh(nx=2, ny=2), order=0)
    assert fes.GetOrder(NodeId(FACE, 3)) == 0
    m = Mesh(unit_square.GenerateMesh(maxh=0.3))
    for space in [HDiv, L2, HDivDiv, HCurlCurl]:
        fes = space(m, order=2)
        assert all(fes.GetOrder(NodeId(FACE, f)) == 2 for f in range(m.nface))
        fes.SetOrder(NodeId(FACE, 5), 3)
        assert [fes.GetOrder(NodeId(FACE, f)) for f in range(m.nface)].count(3) == 1
    m3 = MakeStructured3DMesh(hexes=True, nx=1, ny=2, nz=2)
    fes = SurfaceL2(m3, order=2, definedon=m3.Boundaries("back|bottom"))
    for f in range(m3.nface):
        fes.SetOrder(NodeId(FACE, f), 3)
        fes.GetOrder(NodeId(FACE, f))


def test_order_cap():
    m = MakeStructured2DMesh(nx=2, ny=2)
    for o in [256, 2**15, 2**31-1, 10**9, 1e12]:
        with pytest.raises(Exception):
            H1(m, order=o)
    with pytest.raises(Exception):
        H1(m, order=2, trigorder=300)
    assert H1(m, order=127).globalorder == 127
    fes = H1(m, order=2)
    with pytest.raises(Exception):
        fes.SetOrder(ET.TRIG, 200)
    with pytest.raises(Exception):
        fes.SetOrder(NodeId(EDGE, 0), 1000)
    with pytest.raises(Exception):
        H1(m, dim=10**9)


def test_smoothingblocks_1d():
    fes = H1(Make1DMesh(3), order=3)
    blocks = fes.CreateSmoothingBlocks()
    assert sum(len(b) for b in blocks) == fes.ndof


# ------------------------------------------------------------ mesh bindings

def test_mesh_bindings():
    m = MakeStructured2DMesh(nx=2, ny=2)
    with pytest.raises(Exception):
        m.GetParentFaces(0)
    with pytest.raises(Exception):
        m.LocalHCF()
    with pytest.raises(Exception):
        Make1DMesh(4).SetElementOrder(ElementId(BND, 1), 2)
    with pytest.raises(Exception):
        m.SetElementOrder(ElementId(BND, 1), 2)
    for o in [0, -1, 100, 10**4]:
        with pytest.raises(Exception):
            m.Curve(o)
    assert abs(Integrate(1, m) - 1) < 1e-12
    for o in [-1, 10**7, 2**31-1]:
        with pytest.raises(Exception):
            Integrate(x, Make1DMesh(3), order=o)


def test_interpolate_unknown_operator():
    fes = H1(MakeStructured2DMesh(nx=2, ny=2), order=1)
    with pytest.raises(Exception):
        Interpolate(fes.TestFunction(), fes, operator="bogus")


def test_convert_l2_operator():
    m = MakeStructured2DMesh(nx=2, ny=2)
    fes = FacetFESpace(m, order=1)
    with pytest.raises(Exception):
        fes.ConvertL2Operator(fes)
    h1 = H1(m, order=2)
    with pytest.raises(Exception):
        h1.ConvertL2Operator(L2(MakeStructured2DMesh(nx=2, ny=2), order=2))
    l2 = L2(m, order=2)
    gh1 = GridFunction(h1)
    gh1.Set(x*y)
    gl2 = GridFunction(l2)
    gl2.vec.data = h1.ConvertL2Operator(l2) * gh1.vec
    assert l2err(gl2, x*y, m) < 1e-12


def test_trace_operator_incompatible():
    m = MakeStructured3DMesh(hexes=False, nx=1, ny=1, nz=1)
    fes = L2(m, order=1, dgjumps=True)
    with pytest.raises(Exception):
        fes.TraceOperator(fes, True)
    fes.TraceOperator(FacetFESpace(m, order=1), True)


def test_foreign_regions():
    m = MakeStructured2DMesh(nx=2, ny=2)
    m2 = MakeStructuredSurfaceMesh(nx=2, ny=2)
    with pytest.raises(Exception):
        PlateauFESpace(L2(m2, order=1), [m.BBoundaries(".*")])
    fes = H1(m, order=1)
    with pytest.raises(Exception):
        fes.GetDofs(m2.Boundaries(".*"))
    with pytest.raises(Exception):
        fes.Elements(m2.Boundaries(".*"))
    assert sum(fes.GetDofs(m.Boundaries(".*"))) == 8


def test_from_archive_invalid():
    assert not hasattr(ngsolve.comp, "FromArchiveMesh")
    with pytest.raises(Exception):
        ngsolve.comp.FromArchiveCF("abc", True)
    cf = ngsolve.comp.FromArchiveCF(ngsolve.comp.ToArchive(x*y, True), True)
    assert cf(MakeStructured2DMesh(nx=2, ny=2)(0.5, 0.5)) == pytest.approx(0.25)


def test_pickle_wrapper_spaces():
    m = MakeStructured2DMesh(nx=2, ny=2)
    for fes in [Hidden(H1(m, order=2)), MatrixValued(H1(m, order=2), symmetric=True),
                Reorder(H1(m, order=2)), PlateauFESpace(H1(m, order=2), [m.Boundaries("left")]),
                ngsolve.comp.GlobalInterfaceSpace(m, mapping=x, definedon=m.Boundaries("left"), order=3)]:
        fes2 = pickle.loads(pickle.dumps(fes))
        assert type(fes2) == type(fes)
        assert fes2.ndof == fes.ndof


def test_huge_sizes():
    m = MakeStructured2DMesh(nx=2, ny=2)
    with pytest.raises(Exception):
        L2(m, order=10**4)
    with pytest.raises(Exception):
        GridFunction(H1(m), multidim=10**9)
    with pytest.raises(Exception):
        GridFunction(H1(m), multidim=10**12)
    with pytest.raises(Exception):
        SetHeapSize(10**15)
    assert abs(Integrate(x, m, order=5) - 0.5) < 1e-12


# ------------------------------------------------------------ preconditioner bindings

def test_bddc_class_real_form_complex_space():
    m = MakeStructured2DMesh(nx=4, ny=4)
    fes = H1(m, order=3, dirichlet=".*", complex=True)
    u, v = fes.TnT()
    a = BilinearForm(grad(u)*grad(v)*dx, real=True)
    c = ngsolve.comp.BDDCPreconditioner(a)
    a.Assemble()
    f = LinearForm(1j*v*dx).Assemble()
    gfu = GridFunction(fes)
    gfu.vec.data = solvers.CG(mat=a.mat, pre=c, rhs=f.vec, tol=1e-12, printrates=False)
    ref = GridFunction(fes)
    ref.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
    diff = gfu.vec.CreateVector()
    diff.data = gfu.vec - ref.vec
    assert Norm(diff) < 1e-8 * Norm(ref.vec)


def test_multigrid_preconditioner_old_level():
    m = Mesh(unit_square.GenerateMesh(maxh=0.4))
    fes = H1(m, order=1, dirichlet=".*", autoupdate=True)
    u, v = fes.TnT()
    a = BilinearForm(grad(u)*grad(v)*dx, autoupdate=True).Assemble()
    m.Refine()
    ngsolve.comp.MultiGridPreconditioner(a)


def test_python_preconditioner_creator():
    ngsolve.comp.RegisterPreconditioner("fuzz_comp_pre", lambda mat, freedofs, flags: IdentityMatrix(mat.height))
    m = MakeStructured2DMesh(nx=2, ny=2)
    fes = H1(m, order=1)
    u, v = fes.TnT()
    a = BilinearForm(u*v*dx)
    c = Preconditioner(a, "fuzz_comp_pre")
    x = BaseVector(fes.ndof); x[:] = 1
    with pytest.raises(Exception):
        y = (c * x).Evaluate()
    a.Assemble()
    y = x.CreateVector()
    y.data = c.mat * x
    assert Norm(y) == pytest.approx(Norm(x))
