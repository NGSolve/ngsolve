import pytest
from math import log
from ngsolve import *
from ngsolve.meshes import MakeStructured2DMesh, MakeStructured3DMesh
from netgen.geom2d import unit_square
from netgen.csg import unit_cube

def mesh_of(kind, n):
    if kind == "trig": return Mesh(unit_square.GenerateMesh(maxh=1/n))
    if kind == "quad": return MakeStructured2DMesh(quads=True, nx=n, ny=n)
    if kind == "tet": return Mesh(unit_cube.GenerateMesh(maxh=1/n))
    if kind == "prism": return MakeStructured3DMesh(hexes=False, prism=True, nx=n, ny=n, nz=n)
    return MakeStructured3DMesh(hexes=True, nx=n, ny=n, nz=n)

def stokes_errors(mesh, k, rt, ws):
    # MCS Stokes, sigma = grad u (ws=False) or sigma = 2 eps(u) with weak symmetry
    dim = mesh.dim
    X = (x, y) if dim == 2 else (x, y, z)
    psi = (x*(1-x)*y*(1-y))**2 if dim == 2 else (x*(1-x)*y*(1-y)*z*(1-z))**2
    uex = CF((psi.Diff(y), -psi.Diff(x)) + ((0,) if dim == 3 else ()))
    pex = x-0.5
    gu = CF(tuple(uex[i].Diff(c) for i in range(dim) for c in X), dims=(dim, dim))
    sex = gu+gu.trans if ws else gu
    f = CF(tuple(-sum(sex[i, j].Diff(c) for j, c in enumerate(X)) + pex.Diff(X[i]) for i in range(dim)))

    V = HDiv(mesh, order=k, RT=rt, dirichlet=".*")
    S = V.CreateStressSpace(weaksym=ws)
    sp = [V, S, V.CreateDivergenceSpace()] + ([S.CreateSymmetryMultiplierSpace()] if ws else [])
    fes = FESpace(sp)
    (u, s, p, *e), (v, t, q, *xi) = fes.TnT()
    n = specialcf.normal(dim)
    tg = lambda w: w - (w*n)*n
    B = lambda tau, w: InnerProduct(tau, Grad(w))*dx - (tg(tau*n)*tg(w))*dx(element_boundary=True)
    a = BilinearForm(fes, symmetric=True)
    a += (0.5 if ws else 1)*InnerProduct(s, t)*dx - B(t, u) - B(s, v) - div(u)*q*dx - div(v)*p*dx - 1e-10*p*q*dx
    if ws:
        a += InnerProduct(t, e[0])*dx + InnerProduct(s, xi[0])*dx
    a.Assemble()
    lf = LinearForm(-f*v*dx(bonus_intorder=4)).Assemble()
    gf = GridFunction(fes)
    gf.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="umfpack") * lf.vec
    uh, sh = gf.components[0], gf.components[1]
    assert sqrt(Integrate(div(uh)**2, mesh)) < 1e-8
    sdev = sex - Trace(sex)/dim*Id(dim)
    return (sqrt(Integrate((uh-uex)**2, mesh, order=2*k+6)),
            sqrt(Integrate(InnerProduct(sh-sdev, sh-sdev), mesh, order=2*k+6)))

@pytest.mark.parametrize("kind", ["trig", "quad", "tet", "prism", "hex"])
@pytest.mark.parametrize("k", [1, 2])
@pytest.mark.parametrize("rt", [False, True])
def test_orders(kind, k, rt):
    mesh = mesh_of(kind, 2)
    V = HDiv(mesh, order=k, RT=rt)
    S = V.CreateStressSpace()
    assert type(S) is HCurlDiv
    oi = 2 if (kind == "trig" and rt and k == 1) else k
    assert S.ndof == HCurlDiv(mesh, order=k-1, orderinner=oi).ndof
    if kind != "hex":
        # GG bubbles on simplices replace the RT1 trig orderinner=2
        gg = kind in ("trig", "tet")
        assert V.CreateStressSpace(weaksym=True).ndof == HCurlDiv(mesh, order=k-1, orderinner=k if gg else oi, GGbubbles=gg).ndof
    assert Discontinuous(S).CreateMultiplierSpace().ndof == TangentialFacetFESpace(mesh, order=k-1).ndof
    Sdc = Discontinuous(V.CreateStressSpace())
    assert Sdc.base_space.ndof == S.ndof
    M = Sdc.base_space.CreateSymmetryMultiplierSpace()
    assert type(M) is MatrixValued
    assert M.ndof == MatrixValued(L2(mesh, order=k-1), skewsymmetric=True).ndof
    gf = GridFunction(M)
    gf.vec.SetRandom()
    assert Integrate(Norm(gf+gf.trans), mesh) < 1e-12

def test_explicit_order_and_errors():
    mesh = mesh_of("trig", 2)
    assert HDiv(mesh, order=2).CreateStressSpace(order=2).ndof == HCurlDiv(mesh, order=2).ndof
    with pytest.raises(Exception):
        HDiv(mesh, order=0).CreateStressSpace()

@pytest.mark.parametrize("kind,ns", [("trig", [4, 8]), ("quad", [4, 8]), ("tet", [2, 4]), ("prism", [4, 8]), ("hex", [2, 4])])
@pytest.mark.parametrize("k", [1, 2])
@pytest.mark.parametrize("rt", [False, True])
@pytest.mark.parametrize("ws", [False, True])
def test_stokes_rates(kind, ns, k, rt, ws):
    if kind == "prism" and k == 2:
        ns = [2, 4]  # already asymptotic, n=8 takes a minute
    errs = [stokes_errors(mesh_of(kind, n), k, rt, ws) for n in ns]
    ru = log(errs[0][0]/errs[1][0])/log(ns[1]/ns[0])
    rs = log(errs[0][1]/errs[1][1])/log(ns[1]/ns[0])
    if kind in ("trig", "quad", "prism"):
        assert ru > k+0.5 and rs > k-0.4, (errs, ru, rs)
    else:  # 3D meshes are pre-asymptotic
        assert ru > 1 and rs > 0.2, (errs, ru, rs)

@pytest.mark.parametrize("seed", range(4))
@pytest.mark.parametrize("k", [1, 2, 3])
def test_tet_prism_nt_continuity(seed, k):
    # prism and tet sharing the face z=1: sigma_nt from both sides agrees
    import random
    from netgen.meshing import Mesh as NM, Element3D, MeshPoint, Pnt, FaceDescriptor
    random.seed(seed)
    P = [(0,0,0), (1,0,0), (0,1,0), (0,0,1), (1,0,1), (0,1,1), (0.2,0.3,2)]
    m = NM(dim=3)
    m.Add(FaceDescriptor(surfnr=1, domin=1, bc=1))
    pts = {}
    for i in random.sample(range(7), 7):
        pts[i] = m.Add(MeshPoint(Pnt(*P[i])))
    m.Add(Element3D(1, [pts[i] for i in range(6)]))
    m.Add(Element3D(1, [pts[i] for i in (3, 4, 5, 6)]))
    mesh = Mesh(m)
    g = GridFunction(HDiv(mesh, order=k).CreateStressSpace())
    g.vec.SetRandom()
    eps = 1e-8
    for a, b in [(0.1, 0.1), (0.3, 0.5), (0.6, 0.2)]:
        lo, hi = g(mesh(a, b, 1-eps)), g(mesh(a, b, 1+eps))
        assert abs(lo[2]-hi[2]) < 1e-5 and abs(lo[5]-hi[5]) < 1e-5

@pytest.mark.parametrize("kind", ["tet", "prism", "warped prism"])
@pytest.mark.parametrize("fi", [(0, 1), (1, 2), (2, 3)])
def test_dual_set_reproduces(kind, fi):
    import numpy as np
    if kind == "tet":
        mesh = mesh_of("tet", 2)
    else:
        warp = (lambda x, y, z: (x+0.1*y*z, y+0.05*x*x, z+0.1*x*y)) if kind == "warped prism" else None
        mesh = MakeStructured3DMesh(hexes=False, prism=True, nx=2, ny=2, nz=2, mapping=warp)
    fes = HCurlDiv(mesh, order=fi[0], orderinner=fi[1])
    g, gs = GridFunction(fes), GridFunction(fes)
    g.vec.FV().NumPy()[:] = np.random.RandomState(1).uniform(-1, 1, fes.ndof)
    gs.Set(g, dual=True)
    assert np.abs(gs.vec.FV().NumPy() - g.vec.FV().NumPy()).max() < 1e-9
