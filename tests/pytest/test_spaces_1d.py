"""
H(div), H(curl), normal facet and surface L2 spaces on 1D meshes.

1D de Rham complex: H1 --d/dx--> L2.
HDiv(order=k) is continuous P_{k+1} (vertex values = normal traces), div maps onto L2(order=k);
HCurl(order=k) is discontinuous P_k (1-forms), grad(H1(order=k+1)) lies in it;
NormalFacetFESpace has one dof per vertex, SurfaceL2 one value per boundary point.
"""
import random
from math import log
import numpy as np
import pytest
from netgen.meshing import Mesh as NGMesh, MeshPoint, Pnt, Element0D, Element1D
import ngsolve
from ngsolve import *


def mesh1d(n=6, mapping=lambda x: x + 0.1*x*(1-x), alternate=True):
    """two materials, elements alternately oriented left-to-right and right-to-left"""
    m = NGMesh(dim=1)
    pids = [m.Add(MeshPoint(Pnt(mapping(i/n), 0, 0))) for i in range(n+1)]
    left = m.AddRegion("leftdom", dim=1)
    right = m.AddRegion("rightdom", dim=1)
    bl = m.AddRegion("left", dim=0)
    br = m.AddRegion("right", dim=0)
    for i in range(n):
        pair = [pids[i], pids[i+1]] if (alternate and i % 2) else [pids[i+1], pids[i]]
        m.Add(Element1D(pair, index=left if i < n//2 else right))
    m.Add(Element0D(pids[0], index=bl))
    m.Add(Element0D(pids[n], index=br))
    return Mesh(m)


def poly(deg, seed=1):
    rnd = random.Random(seed)
    return sum(rnd.uniform(-1, 1) * x**j for j in range(deg+1))


def dense(mat):
    return np.array(mat.ToDense())


SPACES = ["HDiv", "HCurl", "NormalFacetFESpace", "SurfaceL2"]


def make(name, mesh, order, **kwargs):
    return getattr(ngsolve, name)(mesh, order=order, **kwargs)


@pytest.mark.parametrize("name", SPACES)
@pytest.mark.parametrize("order", [0, 1, 3])
@pytest.mark.parametrize("definedon", [False, True])
def test_invariants(name, order, definedon):
    mesh = mesh1d()
    kw = {}
    if definedon:
        kw["definedon"] = mesh.Boundaries("left") if name == "SurfaceL2" else mesh.Materials("leftdom")
    fes = make(name, mesh, order, **kw)
    fes._CheckElementFamilies()
    for vb in (VOL, BND):
        fes._CheckDofConsistency(vb)


def test_ndof():
    mesh = mesh1d(n=6)
    for k in range(4):
        assert HDiv(mesh, order=k).ndof == 7 + 6*k
        assert HDiv(mesh, order=k+1, RT=False).ndof == 7 + 6*k
        assert HCurl(mesh, order=k).ndof == 6*(k+1)
        assert HCurl(mesh, order=k+1, type1=True).ndof == 6*(k+1)
        assert NormalFacetFESpace(mesh, order=k).ndof == 7
        assert SurfaceL2(mesh, order=k).ndof == 2


def check_spd(mat, n, rank=None):
    A = dense(mat)
    assert A.shape == (n, n)
    assert np.allclose(A, A.T, atol=1e-13)
    ev = np.linalg.eigvalsh(A)
    assert ev[0] > -1e-12 * max(1, ev[-1])
    if rank is not None:
        assert np.sum(ev > 1e-10 * ev[-1]) == rank


@pytest.mark.parametrize("order", [0, 2])
def test_mass(order):
    mesh = mesh1d()
    fes = HDiv(mesh, order=order)
    u, v = fes.TnT()
    check_spd(BilinearForm(u*v*dx).Assemble().mat, fes.ndof, fes.ndof)
    check_spd(BilinearForm(u.Trace()*v.Trace()*ds).Assemble().mat, fes.ndof, 2)

    fes = HCurl(mesh, order=order)
    u, v = fes.TnT()
    check_spd(BilinearForm(u*v*dx).Assemble().mat, fes.ndof, fes.ndof)

    fes = NormalFacetFESpace(mesh, order=order)
    u, v = fes.TnT()
    check_spd(BilinearForm(u*v*dx(element_boundary=True)).Assemble().mat, fes.ndof, fes.ndof)
    check_spd(BilinearForm(u.Trace()*v.Trace()*ds).Assemble().mat, fes.ndof, 2)

    fes = SurfaceL2(mesh, order=order)
    u, v = fes.TnT()
    check_spd(BilinearForm(u*v*ds).Assemble().mat, fes.ndof, fes.ndof)
    check_spd(BilinearForm(u*v*dx).Assemble().mat, fes.ndof, 0)


@pytest.mark.parametrize("order", [0, 1, 3])
def test_hdiv_reproduce(order):
    mesh = mesh1d()
    for fes, deg in [(HDiv(mesh, order=order), order+1),
                     (HDiv(mesh, order=order+1, RT=False), order+1)]:
        p = poly(deg)
        cf = CF((p,))
        for setter in ("Set", "Interpolate"):
            gf = GridFunction(fes)
            getattr(gf, setter)(cf)
            assert sqrt(Integrate((gf-cf)**2, mesh)) < 1e-12
            assert sqrt(Integrate((div(gf)-p.Diff(x))**2, mesh)) < 1e-11
            assert sqrt(Integrate((gf.Trace()-cf)**2, mesh, BND)) < 1e-12
            for xp in (0, 0.13, 0.5, 0.77, 1):
                assert abs(gf(mesh(xp))[0] - p(mesh(xp))) < 1e-12
        # boundary values only
        gf = GridFunction(fes)
        gf.Set(cf, definedon=mesh.Boundaries(".*"))
        assert sqrt(Integrate((gf.Trace()-cf)**2, mesh, BND)) < 1e-12
        assert gf.vec.FV().NumPy()[0] == pytest.approx(p(mesh(0.)))


@pytest.mark.parametrize("order", [0, 1, 3])
def test_hdiv_dual_commutes(order):
    # div of the interpolant = L2 projection of div onto L2(order)
    mesh = mesh1d()
    fes = HDiv(mesh, order=order)
    l2 = L2(mesh, order=order)
    cf = CF((sin(3*x)+x**(order+3),))
    gf = GridFunction(fes)
    gf.Interpolate(cf, bonus_intorder=12)
    gl2 = GridFunction(l2)
    gl2.Set(cf[0].Diff(x), bonus_intorder=12)
    assert sqrt(Integrate((div(gf)-gl2)**2, mesh, order=2*order+2)) < 1e-12
    # vertex values are interpolated
    for vt in mesh.vertices:
        assert abs(gf.vec[vt.nr] - cf[0](mesh(vt.point[0]))) < 1e-12


@pytest.mark.parametrize("order", [0, 1, 3])
def test_hcurl_reproduce(order):
    mesh = mesh1d()
    for fes, deg in [(HCurl(mesh, order=order), order),
                     (HCurl(mesh, order=order+1, type1=True), order)]:
        p = poly(deg, seed=2)
        cf = CF((p,))
        for setter in ("Set", "Interpolate"):
            gf = GridFunction(fes)
            getattr(gf, setter)(cf)
            assert sqrt(Integrate((gf-cf)**2, mesh)) < 1e-12
            for xp in (0.13, 0.5, 0.77):
                assert abs(gf(mesh(xp))[0] - p(mesh(xp))) < 1e-12


@pytest.mark.parametrize("order", [0, 1, 3])
def test_grad_h1_in_hcurl(order):
    mesh = mesh1d()
    h1 = H1(mesh, order=order+1)
    gh = GridFunction(h1)
    gh.vec.FV().NumPy()[:] = np.random.default_rng(3).uniform(-1, 1, h1.ndof)
    fes = HCurl(mesh, order=order)
    for setter in ("Set", "Interpolate"):
        gf = GridFunction(fes)
        getattr(gf, setter)(grad(gh))
        assert sqrt(Integrate((gf-grad(gh))**2, mesh)) < 1e-12


def test_normalfacet_reproduce():
    mesh = mesh1d()
    fes = NormalFacetFESpace(mesh, order=2)
    cf = CF((1+x**3,))
    n = specialcf.normal(1)
    for setter in ("Set", "Interpolate"):
        gf = GridFunction(fes)
        if setter == "Set":
            gf.Set(cf, dual=True)
        else:
            gf.Interpolate(cf)
        assert sqrt(Integrate((gf-cf)**2*dx(element_boundary=True), mesh)) < 1e-12
        assert sqrt(Integrate((gf.Operator("normalcomponent")-cf*n)**2*dx(element_boundary=True), mesh)) < 1e-12
        for v in mesh.vertices:
            assert abs(gf.vec[v.nr] - cf[0](mesh(v.point[0]))) < 1e-12
    gf = GridFunction(fes)
    gf.Set(cf, definedon=mesh.Boundaries(".*"))
    assert sqrt(Integrate((gf-cf)**2, mesh, BND)) < 1e-12


def test_surfacel2_reproduce():
    mesh = mesh1d()
    fes = SurfaceL2(mesh, order=2)
    cf = 2+x
    for setter in ("Set", "Interpolate"):
        gf = GridFunction(fes)
        getattr(gf, setter)(cf, definedon=mesh.Boundaries(".*"))
        assert sqrt(Integrate((gf-cf)**2, mesh, BND)) < 1e-12
        vals = Integrate(gf, mesh, BND, region_wise=True)
        assert abs(vals[0] - 2) < 1e-12 and abs(vals[1] - 3) < 1e-12


def mixed_poisson(mesh, k, f, RT=True):
    order = k if RT else k+1
    V = HDiv(mesh, order=order, RT=RT)
    Q = L2(mesh, order=k)
    sigma, u = (V*Q).TrialFunction()
    tau, v = (V*Q).TestFunction()
    a = BilinearForm((sigma*tau + div(sigma)*v + div(tau)*u)*dx).Assemble()
    rhs = LinearForm(-f*v*dx).Assemble()
    gf = GridFunction(V*Q)
    gf.vec.data = a.mat.Inverse(inverse="umfpack") * rhs.vec
    return gf.components


@pytest.mark.parametrize("k", [0, 1, 2])
def test_mixed_poisson_convergence(k):
    uex = sin(pi*x)
    f = pi**2*sin(pi*x)
    errs = []
    for n in (4, 8, 16):
        mesh = mesh1d(n=n)
        sigma, u = mixed_poisson(mesh, k, f)
        errs.append(sqrt(Integrate((u-uex)**2, mesh)))
    rates = [log(errs[i]/errs[i+1])/log(2) for i in range(2)]
    assert min(rates) > k+1-0.15


@pytest.mark.parametrize("k", [0, 1, 2])
def test_mixed_poisson_exact(k):
    # f in P_k: sigma_h = sigma exactly (sigma_h - sigma is constant with mean 0),
    # u_h = L2 projection of u onto L2(order=k)
    rnd = random.Random(5)
    a = [rnd.uniform(-1, 1) for j in range(k+1)]
    f = sum(a[j]*x**j for j in range(k+1))
    c1 = sum(a[j]/((j+1)*(j+2)) for j in range(k+1))
    uex = -sum(a[j]*x**(j+2)/((j+1)*(j+2)) for j in range(k+1)) + c1*x
    mesh = mesh1d(n=5)
    sigma, u = mixed_poisson(mesh, k, f)
    assert sqrt(Integrate((sigma-uex.Diff(x))**2, mesh)) < 1e-12
    pu = GridFunction(L2(mesh, order=k))
    pu.Set(uex, bonus_intorder=4)
    assert sqrt(Integrate((u-pu)**2, mesh)) < 1e-12


@pytest.mark.parametrize("k", [0, 1, 3])
def test_div_onto_l2(k):
    mesh = mesh1d()
    for V in (HDiv(mesh, order=k), HDiv(mesh, order=k+1, RT=False)):
        Q = L2(mesh, order=k)
        B = BilinearForm(div(V.TrialFunction())*Q.TestFunction()*dx).Assemble()
        Bd = dense(B.mat)
        assert Bd.shape == (Q.ndof, V.ndof)
        assert np.linalg.matrix_rank(Bd) == Q.ndof
    # kernel of div: the global constants
    V = HDiv(mesh, order=k)
    assert V.ndof - Q.ndof == 1
    assert HDiv(mesh, order=k).CreateDivergenceSpace().ndof == Q.ndof


@pytest.mark.parametrize("k", [1, 2])
def test_normalfacet_hybrid_poisson(k):
    # primal hybrid method: u in L2, lam = normal flux in NormalFacetFESpace
    uex = sin(pi*x)
    f = pi**2*sin(pi*x)
    errs = []
    for nel in (4, 8, 16):
        mesh = mesh1d(n=nel)
        X = L2(mesh, order=k) * NormalFacetFESpace(mesh)
        (u, lam), (v, mu) = X.TnT()
        n = specialcf.normal(1)
        a = BilinearForm(grad(u)*grad(v)*dx - (lam*n)*v*dx(element_boundary=True)
                         - (mu*n)*u*dx(element_boundary=True)).Assemble()
        rhs = LinearForm(f*v*dx).Assemble()
        gf = GridFunction(X)
        gf.vec.data = a.mat.Inverse(inverse="umfpack") * rhs.vec
        errs.append(sqrt(Integrate((gf.components[0]-uex)**2, mesh)))
        # the multiplier is the flux u' at the vertices
        flux = [abs(gf.components[1].vec[vt.nr] - pi*cos(pi*vt.point[0])) for vt in mesh.vertices]
        assert max(flux) < 2 * (1/nel)**k
    rates = [log(errs[i]/errs[i+1])/log(2) for i in range(2)]
    assert min(rates) > k+1-0.15


def test_rejected_spaces():
    mesh = mesh1d()
    with pytest.raises(Exception, match="TangentialFacetFESpace is not defined on 1D meshes"):
        TangentialFacetFESpace(mesh, order=1)
    with pytest.raises(Exception, match="VectorSurfaceL2 is not defined on 1D meshes"):
        VectorSurfaceL2(mesh, order=1)


def test_context_gate():
    mesh = mesh1d()
    u, v = HCurl(mesh, order=1).TnT()
    with pytest.raises(Exception):
        BilinearForm(u*v*ds).Assemble()
    with pytest.raises(Exception):
        curl(u)
    lam, mu = NormalFacetFESpace(mesh).TnT()
    with pytest.raises(Exception, match="inside elements"):
        BilinearForm(lam*mu*dx).Assemble()
    gf = GridFunction(NormalFacetFESpace(mesh))
    with pytest.raises(Exception):
        Integrate(gf*gf, mesh)
    # accepted
    s, t = HDiv(mesh, order=1).TnT()
    BilinearForm(s*t*dx + div(s)*div(t)*dx + s.Trace()*t.Trace()*ds + s.Operator("normalcomponent")*t.Operator("normalcomponent")*dx(element_boundary=True)).Assemble()
    BilinearForm(InnerProduct(grad(s), grad(t))*dx).Assemble()
    BilinearForm(lam*mu*dx(element_boundary=True) + lam.Trace()*mu.Trace()*ds).Assemble()
