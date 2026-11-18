"""
Illegal constructs are rejected when they are built or bound, with a clean
exception instead of a crash; legal code keeps working.
"""
import gc
import pickle
import pytest
import numpy as np
from netgen.geom2d import unit_square, SplineGeometry
from ngsolve import *
from ngsolve.fem import BFI, LFI, H1FE, MixedFE
import ngsolve.comp


@pytest.fixture(scope="module")
def m():
    return Mesh(unit_square.GenerateMesh(maxh=0.3))


def fresh_mesh(maxh=0.3):
    return Mesh(unit_square.GenerateMesh(maxh=maxh))


# ---------------------------------------------------------------- fem

def test_placeholder_cycle(m):
    p = PlaceholderCF(x)
    with pytest.raises(Exception, match="placeholder itself"):
        p.Set(p + 1)
    p.Set(y + 1)
    assert Integrate(p, m) == pytest.approx(1.5)


def test_bfi_too_few_coefs(m):
    with pytest.raises(Exception, match="coefficient components"):
        BFI("elasticity", coef=1, dim=2)


def test_bfi_vector_coef(m):
    with pytest.raises(Exception, match="3 components"):
        BFI("mass", coef=CF((1, 2, 3)))


def test_named_bfi_valid(m):
    fes = H1(m, order=1)
    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += BFI("mass", coef=2)
    a.Assemble()
    b = BilinearForm(2 * u * v * dx).Assemble()
    w = a.mat.CreateColVector()
    gf = GridFunction(fes)
    gf.Set(x)
    w.data = a.mat * gf.vec - b.mat * gf.vec
    assert Norm(w) < 1e-12


def test_named_lfi_wrong_space(m):
    fes = HCurl(m, order=1)
    f = LinearForm(fes)
    with pytest.raises(Exception, match="does not match"):
        f += LFI("source", coef=1)


def test_named_lfi_valid(m):
    fes = H1(m, order=1)
    f = LinearForm(fes)
    f += LFI("source", coef=1)
    f.Assemble()
    assert sum(f.vec) == pytest.approx(1)


def test_named_bfi_space_dimension(m):
    a = BilinearForm(H1(m, order=1))
    with pytest.raises(Exception, match="dimension"):
        a += BFI("elasticity", coef=[1, 0.3])
    a = BilinearForm(H1(m, order=1, dim=2))
    a += BFI("elasticity", coef=[1, 0.3])
    a.Assemble()
    assert a.mat.height == m.nv


@pytest.mark.parametrize("dim,order", [(2, 1), (2, 3), (3, 1), (3, 2)])
def test_named_bfi_elasticity_values(dim, order):
    if dim == 2:
        mesh = fresh_mesh()
    else:
        from netgen.csg import unit_cube
        mesh = Mesh(unit_cube.GenerateMesh(maxh=0.5))
    E, nu = 2.5, 0.3
    mu, lam = E / 2 / (1 + nu), E * nu / ((1 + nu) * (1 - 2 * nu))
    fes = H1(mesh, order=order, dim=dim)
    a = BilinearForm(fes)
    a += BFI("elasticity", coef=[E, nu])
    a.Assemble()
    u, v = fes.TnT()
    b = BilinearForm(2 * mu * InnerProduct(Sym(Grad(u)), Sym(Grad(v))) * dx
                     + lam * Trace(Grad(u)) * Trace(Grad(v)) * dx).Assemble()
    da, db = np.array(a.mat.ToDense()), np.array(b.mat.ToDense())
    assert np.max(np.abs(da - db)) < 1e-10 * np.max(np.abs(db))


def test_mixedfe_mismatch(m):
    with pytest.raises(Exception, match="MixedFE"):
        MixedFE(H1FE(ET.TRIG, 1), H1FE(ET.TET, 2))


def test_mixedfe_keeps_elements(m):
    fe = MixedFE(H1FE(ET.TRIG, 1), H1FE(ET.TRIG, 2))
    gc.collect()
    u, v = H1(m, order=2).TnT()
    mat = SymbolicBFI(u * v).CalcElementMatrix(fe, m.GetTrafo(ElementId(VOL, 0)))
    assert mat.h == 6 and mat.w == 3


def test_pickle_data_tamper(m):
    c = x * y
    d = c.data
    d["data"] = d["data"][:len(d["data"]) // 2]
    with pytest.raises(Exception, match="corrupted"):
        CoefficientFunction(d)
    c2 = CoefficientFunction(c.data)
    assert c2(m(0.2, 0.3)) == pytest.approx(0.06)
    c3 = pickle.loads(pickle.dumps(c))
    assert c3(m(0.2, 0.3)) == pytest.approx(0.06)


def test_meshpoint_forged_nr(m):
    pts = m(np.array([0.2, 0.3]), np.array([0.2, 0.3]))
    assert np.allclose(x(pts).flatten(), [0.2, 0.3])
    pts2 = pts.copy()
    pts2["nr"] = 10**6
    with pytest.raises(IndexError):
        x(pts2)
    pts3 = pts.copy()
    pts3["VorB"] = 77
    with pytest.raises(Exception, match="VorB"):
        x(pts3)
    pts4 = pts.copy()
    pts4["meshptr"] = 4096
    with pytest.raises(Exception, match="mesh"):
        x(pts4)


def test_meshpoint_dead_mesh():
    mesh = fresh_mesh()
    pts = mesh(np.array([0.2, 0.3]), np.array([0.2, 0.3])).copy()
    del mesh
    gc.collect()
    with pytest.raises(Exception, match="mesh"):
        x(pts)


def test_trafo_mip_keeps_mesh():
    def f():
        mesh = fresh_mesh()
        return mesh.GetTrafo(ElementId(VOL, 0))(0.2, 0.2)
    mip = f()
    gc.collect()
    junk = [fresh_mesh(0.5) for i in range(3)]
    val = x(mip)
    assert 0 <= val <= 1


def test_trafo_mip_views_keep_mip(m):
    trafo = m.GetTrafo(ElementId(VOL, 0))
    ref = trafo(0.2, 0.3)
    p, jac = tuple(ref.point), np.array(ref.jacobi)
    assert np.linalg.norm(jac) > 0
    del ref
    gc.collect()
    assert tuple(trafo(0.2, 0.3).point) == pytest.approx(p)
    assert np.array(trafo(0.2, 0.3).jacobi) == pytest.approx(jac)
    view = trafo(0.2, 0.3).point
    junk = [trafo(0.7, 0.1) for i in range(10)]
    assert tuple(view) == pytest.approx(p)
    ir = IntegrationRule([(0.2, 0.3)], [1.])
    assert tuple(trafo(ir[0]).point) == pytest.approx(p)


def test_standalone_trafo_gf(m):
    from ngsolve.fem import ElementTransformation
    t = ElementTransformation(ET.TRIG, [(0, 0), (1, 0), (0, 1)])
    mip = t(0.2, 0.2)
    gf = GridFunction(H1(m, order=2))
    gf.Set(x)
    # a point of no mesh: point search in the mesh of the GridFunction
    assert gf(mip) == pytest.approx(mip.point[0])


def test_meshpoint_keeps_mesh():
    mp = fresh_mesh()(0.3, 0.4)
    gc.collect()
    junk = [fresh_mesh(0.05) for i in range(3)]
    assert x(mp) == pytest.approx(0.3)


def test_domainwise_empty():
    with pytest.raises(Exception, match="empty"):
        IfPos(x, [], 1)
    with pytest.raises(Exception, match="empty"):
        CF([])


def test_domainwise_none_replace(m):
    c = CF([x, None])
    assert c(m(0.2, 0.2)) == pytest.approx(0.2)
    c2 = c.Replace({x: y})
    assert c2(m(0.2, 0.3)) == pytest.approx(0.3)
    assert c.Diff(x)(m(0.2, 0.3)) == pytest.approx(1)
    assert CF([None])(m(0.2, 0.3)) == 0


def test_voxel_trafo_dim(m):
    vals = np.arange(9.).reshape(3, 3)
    with pytest.raises(Exception, match="trafocf"):
        VoxelCoefficient((0, 0), (1, 1), vals, trafocf=CF((x, y, x, y, x, y, x, y)))
    c = VoxelCoefficient((0, 0), (1, 1), vals, trafocf=CF((y, x)))
    assert Integrate(c, m) == pytest.approx(4)


def test_voxel_empty_and_shape(m):
    with pytest.raises(Exception, match="empty"):
        VoxelCoefficient((0, 0), (1, 1), np.ones((0, 0)))
    with pytest.raises(Exception, match="entries"):
        VoxelCoefficient((0, 0), (1, 1, 1), np.ones((2, 2)))
    with pytest.raises(Exception, match="dimensional"):
        VoxelCoefficient((0, 0), (1, 1), np.ones((2, 2, 2)))
    c = VoxelCoefficient((0, 0, 0), (1, 1, 1), np.ones((2, 2, 2)))
    with pytest.raises(Exception, match="dimension"):
        Integrate(c, m)


def test_voxel_nearest_upper_end(m):
    # nearest-voxel lookup at the upper end of a box ending at 0
    c = VoxelCoefficient((-1, -1), (0, 0), np.arange(4.).reshape(2, 2), linear=False)
    mp = m(0.0, 0.0)
    assert c(mp) == pytest.approx(3)


# ---------------------------------------------------------------- comp

def test_order_negative(m):
    for sp in [H1, L2, HCurl, HDiv, FacetFESpace, HCurlDiv, HDivDiv, HCurlCurl]:
        with pytest.raises(Exception, match="order"):
            sp(m, order=-1)
    assert L2(m, order=0).ndof == m.ne


def test_mg_refined_before():
    # coarse matrices come from the Galerkin projection
    mesh = fresh_mesh()
    mesh.Refine()
    mesh.Refine()
    fes = H1(mesh, order=1, dirichlet=".*")
    u, v = fes.TnT()
    a = BilinearForm(grad(u) * grad(v) * dx)
    c = Preconditioner(a, "multigrid")
    a.Assemble()
    f = LinearForm(v * dx).Assemble()
    gfu = GridFunction(fes)
    inv = CGSolver(a.mat, c.mat, tol=1e-10, maxiter=40)
    gfu.vec.data = inv * f.vec
    ref = GridFunction(fes)
    ref.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec
    assert inv.iterations < 25
    assert Integrate((gfu - ref) ** 2, mesh) < 1e-16


def test_mg_nonassemble(m):
    fes = H1(m, order=1)
    u, v = fes.TnT()
    a = BilinearForm(fes, nonassemble=True)
    a += grad(u) * grad(v) * dx + u * v * dx
    with pytest.raises(Exception, match="nonassemble"):
        Preconditioner(a, "multigrid")


def test_mg_valid():
    mesh = fresh_mesh()
    fes = H1(mesh, order=1)
    u, v = fes.TnT()
    a = BilinearForm(grad(u) * grad(v) * dx + u * v * dx)
    c = Preconditioner(a, "multigrid")
    a.Assemble()
    mesh.Refine()
    fes.Update()
    a.Assemble()
    x_ = a.mat.CreateColVector()
    x_[:] = 1
    y_ = (c.mat * x_).Evaluate()
    assert Norm(y_) > 0


def test_prec_before_assemble(m):
    fes = H1(m, order=1)
    u, v = fes.TnT()
    a = BilinearForm(grad(u) * grad(v) * dx + u * v * dx)
    c = Preconditioner(a, "local")
    with pytest.raises(Exception, match="not ready"):
        c.Test()
    with pytest.raises(Exception, match="not ready"):
        c.mat


def test_convert_operator_two_meshes(m):
    m2 = fresh_mesh(0.2)
    with pytest.raises(Exception, match="same mesh"):
        ConvertOperator(H1(m, order=1), H1(m2, order=1))


def test_harmonic_prolongation_other_space():
    mesh = fresh_mesh()
    X = H1(mesh, order=1)
    Y = H1(mesh, order=3)
    u, v = Y.TnT()
    a = BilinearForm(grad(u) * grad(v) * dx)
    with pytest.raises(Exception, match="this space"):
        X.SetHarmonicProlongation(a)


def test_flux_no_integrators(m):
    fes = H1(m, order=2)
    gf = GridFunction(fes)
    a = BilinearForm(fes)
    with pytest.raises(Exception, match="no integrators"):
        a.Flux(gf)


def test_set_preconditioner_none(m):
    fes = H1(m, order=1)
    u, v = fes.TnT()
    a = BilinearForm(u * v * dx)
    a.UnsetPreconditioner(None)
    with pytest.raises(Exception, match="None"):
        a.SetPreconditioner(None)
    a.Assemble()


def test_elements_of_temp_mesh():
    it = fresh_mesh().Elements(VOL)
    gc.collect()
    junk = [fresh_mesh(0.05) for i in range(3)]
    verts = [e.vertices for e in it]
    nv = max(v.nr for vs in verts for v in vs) + 1
    assert nv == fresh_mesh().nv


def test_vertices_of_temp_mesh():
    vs = fresh_mesh().vertices
    gc.collect()
    junk = [fresh_mesh(0.05) for i in range(3)]
    assert vs[0].point == pytest.approx((0, 0))


def test_stale_space_getdofnrs():
    mesh = fresh_mesh()
    fes = H1(mesh, order=2, autoupdate=False)
    mesh.Refine()
    with pytest.raises(Exception, match="not up to date"):
        fes.GetDofNrs(ElementId(VOL, mesh.ne - 1))
    fes.Update()
    assert len(fes.GetDofNrs(ElementId(VOL, mesh.ne - 1))) == 6


def test_tp_helpers_on_ordinary_spaces(m):
    fes = H1(m, order=2)
    gf = GridFunction(fes)
    with pytest.raises(Exception, match="two spaces"):
        ngsolve.comp.TensorProductFESpace([fes])
    with pytest.raises(Exception, match="TensorProductFESpace"):
        ngsolve.comp.Prolongate(gf, gf)
    with pytest.raises(Exception, match="TensorProductFESpace"):
        ngsolve.comp.Transfer2StdMesh(gf, gf)
    with pytest.raises(Exception, match="TensorProductFESpace"):
        ngsolve.comp.TensorProductIntegrate(gf, [0.5, 0.5], CF(1))


# ---------------------------------------------------------------- spaces: GetFE and GetDofNrs agree

from test_element_family import CASES as FAMILY_CASES, SPACES as FAMILY_SPACES, get_mesh as family_mesh


@pytest.mark.parametrize("space,meshname,restricted", FAMILY_CASES,
                         ids=[f"{s}-{mn}-{'definedon' if r else 'all'}" for s, mn, r in FAMILY_CASES])
def test_dof_consistency(space, meshname, restricted):
    mesh = family_mesh(meshname)
    try:
        fes = FAMILY_SPACES[space][0](mesh, restricted)
    except Exception as e:
        pytest.skip(f"construction rejected: {e}")
    for vb in [VOL, BND, BBND, BBBND][:mesh.dim + 1]:
        try:
            fes._CheckDofConsistency(vb)
        except Exception as e:
            # elements a space does not provide at all are fine, inconsistent ones are not
            assert "inconsistent" not in str(e), f"{vb}: {e}"


def test_hdivdivfacet_dofs_distinct():
    mesh = family_mesh("tet3d")
    fes = HDivDivFacetSpace(mesh, order=2)
    used = set()
    for el in fes.Elements(VOL):
        used.update(el.dofs)
    assert len(used) == fes.ndof


# ---------------------------------------------------------------- out-of-bounds accesses found with range checks

def test_todense_block_matrix(m):
    fes = H1(m, order=1, dim=2)
    u, v = fes.TnT()
    a = BilinearForm(InnerProduct(u, v) * dx).Assemble()
    d = np.array(a.mat.ToDense())
    assert d.shape == (2 * m.nv, 2 * m.nv)
    gf = GridFunction(fes)
    gf.Set(CF((x, y)))
    w = a.mat.CreateColVector()
    w.data = a.mat * gf.vec
    assert np.allclose(d @ gf.vec.FV().NumPy(), w.FV().NumPy())


def test_named_bfi_calcelementmatrix_dim(m):
    trafo = m.GetTrafo(ElementId(VOL, 0))
    fe = H1FE(ET.TRIG, 1)
    mat = BFI("elasticity", coef=[1, 0.3]).CalcElementMatrix(fe, trafo)
    ref = BFI("elasticity", coef=[1, 0.3], dim=2).CalcElementMatrix(fe, trafo)
    assert mat.h == 6 and mat.w == 6
    assert np.allclose(np.array(mat), np.array(ref))
    vec = LFI("source", coef=1).CalcElementVector(fe, trafo)
    assert len(vec) == 3


def test_symbolic_bfi_wrong_element(m):
    u, v = HCurl(m, order=2).TnT()
    with pytest.raises(Exception, match="expected"):
        SymbolicBFI(curl(u) * curl(v)).CalcElementMatrix(H1FE(ET.TRIG, 2), m.GetTrafo(ElementId(VOL, 0)))


def test_einsum_scalar_inputs(m):
    from ngsolve.fem import Einsum
    c = Einsum("i,i->", x, y)
    assert c(m(0.2, 0.3)) == pytest.approx(0.06)


def test_cf_empty_slice(m):
    d = CF((1, 2, 3, 4))[3:1]
    assert d.dim == 0
    assert d(m(0.2, 0.2)) == ()


def test_eig_needs_square_matrix(m):
    with pytest.raises(Exception, match="square"):
        CF((1, 2, 3)).Eig()
    e = CF((2, 0, 0, 3), dims=(2, 2)).Eig()
    assert e(m(0.2, 0.2))[4:] == pytest.approx((2, 3))


def test_definedonelements_short(m):
    fes = H1(m, order=1)
    u, v = fes.TnT()
    a = BilinearForm(fes)
    with pytest.raises(Exception, match="definedonelements"):
        a += BFI("laplace", coef=1, definedonelem=BitArray(3))
    with pytest.raises(Exception, match="definedonelements"):
        BilinearForm(u * v * dx(definedonelements=BitArray(3)))
    mesh = fresh_mesh()
    els = BitArray(mesh.ne)
    els.Set()
    mesh.Refine()
    with pytest.raises(Exception, match="definedonelements"):
        Integrate(CF(1) * dx(definedonelements=els), mesh)


def test_hcurlamg_needs_lowest_order_hcurl(m):
    fes = H1(m, order=1)
    u, v = fes.TnT()
    a = BilinearForm(grad(u) * grad(v) * dx + u * v * dx)
    c = Preconditioner(a, "hcurlamg")
    with pytest.raises(Exception, match="hcurlamg"):
        a.Assemble()


def test_stale_spaces_rejected():
    mesh = fresh_mesh()
    f1 = H1(mesh, order=2, autoupdate=False)
    u, v = f1.TnT()
    mesh.Refine()
    with pytest.raises(Exception, match="not up to date"):
        BilinearForm(u * v * dx).Assemble()

    mesh = fresh_mesh()
    f1 = H1(mesh, order=2, autoupdate=False)
    g1 = GridFunction(f1, autoupdate=False)
    mesh.Refine()
    f1.Update()
    with pytest.raises(Exception, match="not up to date"):
        g1.Set(x)

    mesh = fresh_mesh()
    f1 = H1(mesh, order=2, autoupdate=False)
    f2 = L2(mesh, order=1, autoupdate=False)
    X = f1 * f2
    mesh.Refine()
    f1.Update()
    (u, p), (v, q) = X.TnT()
    with pytest.raises(Exception, match="not up to date"):
        BilinearForm(u * v * dx + p * q * dx).Assemble()

    mesh = fresh_mesh()
    fd = VectorH1(mesh, order=1, autoupdate=False)
    d = GridFunction(fd, autoupdate=False)
    d.Set(CF((x, y)))
    mesh.SetDeformation(d)
    mesh.Refine()
    with pytest.raises(Exception, match="not up to date"):
        Integrate(CF(1), mesh)


def test_convert_operator_range_dofs_short(m):
    with pytest.raises(Exception, match="range_dofs"):
        ConvertOperator(H1(m, order=1), L2(m, order=1), range_dofs=BitArray(3))


def test_interpolate_dimension(m):
    with pytest.raises(Exception, match="dimension"):
        Interpolate(CF((x, y)), H1(m, order=2))
    I = Interpolate(CF((x, y)), VectorH1(m, order=2))
    assert Integrate(I, m) == pytest.approx((0.5, 0.5))


def test_contact_minion_must_be_boundary(m):
    with pytest.raises(Exception, match="boundary region"):
        ContactBoundary(m.Materials(".*"), m.Materials(".*"))


def test_regions_of_other_mesh():
    g = SplineGeometry()
    g.AddRectangle((0, 0), (1, 1), bcs=["a", "b", "c", "d"])
    g.AddCircle((0.5, 0.5), 0.2, leftdomain=2, rightdomain=1, bc="circ")
    mb = Mesh(g.GenerateMesh(maxh=0.2))
    m1 = fresh_mesh()
    with pytest.raises(Exception, match="another mesh"):
        H1(mb, order=1, dirichlet=m1.Boundaries("left"))
    with pytest.raises(Exception, match="another mesh"):
        H1(mb, order=1, definedon=m1.Materials("default"))
    with pytest.raises(Exception, match="BND"):
        H1(mb, order=1, dirichlet=mb.Materials(".*"))
    with pytest.raises(Exception, match="another mesh"):
        Integrate(CF(1), mb, definedon=m1.Materials("default"))
    gf = GridFunction(H1(mb, order=1))
    with pytest.raises(Exception, match="another mesh"):
        gf.Set(1, definedon=m1.Boundaries("left|right|top|bottom"))


def test_space_dim_positive(m):
    with pytest.raises(Exception, match="dim"):
        H1(m, order=2, dim=0)
    with pytest.raises(Exception, match="dim"):
        H1(m, order=2, dim=-2)


def test_domainwise_dimension(m):
    # shapes () and (1,) are compatible, only the total dimension must agree
    cf = m.MaterialCF({".*": CF((x,))}, default=CF(1))
    assert Integrate(cf, m) == pytest.approx(Integrate(x, m))
    with pytest.raises(Exception, match="same dimension"):
        CoefficientFunction([CF((x, y)), CF(1)])


def test_freedofs_of_finer_level():
    # hierarchical dof numbering: fine-level freedofs are used for coarse matrices
    from ngsolve.meshes import MakeStructured2DMesh
    mesh = MakeStructured2DMesh(quads=False, nx=2, ny=2)
    fes = H1(mesh, order=1, dirichlet="left")
    u, v = fes.TnT()
    a0 = BilinearForm(grad(u)*grad(v)*dx + u*v*dx).Assemble()
    coarsemat = a0.mat
    mesh.Refine()
    fes.Update()
    assert fes.FreeDofs().__len__() > coarsemat.height
    coarsemat.CreateSmoother(fes.FreeDofs())
    inv = coarsemat.Inverse(fes.FreeDofs())
    w = coarsemat.CreateColVector()
    w[:] = 1
    assert (inv * w).Evaluate().Norm() > 0
    short = BitArray(coarsemat.height - 1)
    with pytest.raises(Exception, match="at least"):
        coarsemat.CreateSmoother(short)
    with pytest.raises(Exception, match="at least"):
        coarsemat.Inverse(short)


def test_linearization_lowest_order_space_with_dim():
    # low-order restriction of the linearization point keeps the entry size
    from ngsolve.meshes import MakeStructured2DMesh
    mesh = MakeStructured2DMesh(quads=False, nx=4, ny=4)
    fes = H1(mesh, order=2, dirichlet="left", dim=2)
    u = fes.TrialFunction()
    a = BilinearForm(fes, symmetric=True)
    a += Variation((InnerProduct(Grad(u), Grad(u)) + InnerProduct(u, u)**2) * dx)
    gfu = GridFunction(fes)
    gfu.vec[:] = 0.1
    a.Assemble()
    pre = Preconditioner(a, "multigrid")
    a.AssembleLinearization(gfu.vec)
    w = gfu.vec.CreateVector()
    w.data = pre.mat * gfu.vec
    assert w.Norm() > 0


def test_other_bnd_with_trial_function():
    # boundary values of Other() may depend on the trial-function (DG apply)
    from ngsolve.meshes import MakeStructured2DMesh
    mesh = MakeStructured2DMesh(quads=False, nx=3, ny=3)
    fes = L2(mesh, order=1)**2
    U, V = fes.TnT()
    h, hu = U
    ubnd = CF((mesh.BoundaryCF({"left": h}, default=0), mesh.BoundaryCF({"left": -hu}, default=0)))
    n = specialcf.normal(2)
    a = BilinearForm(fes, nonassemble=True)
    a += InnerProduct(U - U.Other(bnd=ubnd), V) * dx(element_boundary=True)
    gfu = GridFunction(fes)
    gfu.Set(CF((1+x, y)))
    res = gfu.vec.CreateVector()
    a.Apply(gfu.vec, res)
    assert res.Norm() > 0
    b = BilinearForm(fes, nonassemble=True)
    with pytest.raises(Exception, match="trial-functions only"):
        b += InnerProduct(U - U.Other(bnd=CF((V[0], 0))), V) * dx(element_boundary=True)
        b.Apply(gfu.vec, res)
    c = BilinearForm(fes, nonassemble=True)
    with pytest.raises(Exception, match="trial-functions only"):
        c += InnerProduct(U - U.Other(bnd=U.Other()), V) * dx(element_boundary=True)
        c.Apply(gfu.vec, res)
    # a trial-function of an incompatible space
    other = HDiv(mesh, order=1)
    d = BilinearForm(fes, nonassemble=True)
    with pytest.raises(Exception):
        d += InnerProduct(U - U.Other(bnd=other.TrialFunction()), V) * dx(element_boundary=True)
        d.Apply(gfu.vec, res)


def test_webgui_deformation_without_trace():
    # L2 functions have no evaluator on the mesh edges, the edges fall back
    webgui = pytest.importorskip("ngsolve.webgui")
    from ngsolve.meshes import MakeStructured2DMesh, MakeStructured3DMesh
    for mesh in [MakeStructured2DMesh(quads=False, nx=2, ny=2), MakeStructured3DMesh(hexes=False, nx=1, ny=1, nz=1)]:
        gf = GridFunction(L2(mesh, order=1))
        gf.Set(x)
        regions = {vb: mesh.Region(vb, ".*") for vb in [VOL, BND, BBND, BBBND]}
        d = webgui.BuildRenderData(mesh, gf, deformation=CF((0, 0, 0.02*gf)), regions=regions)
        assert len(d["edges"]) > 0
