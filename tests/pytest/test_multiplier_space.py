import pytest
import numpy as np
from ngsolve import *
from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from ngsolve.meshes import MakeStructured2DMesh, MakeStructured3DMesh

meshes = {
    "trig": lambda: Mesh(unit_square.GenerateMesh(maxh=0.4)),
    "quad": lambda: MakeStructured2DMesh(quads=True, nx=3, ny=3),
    "tet": lambda: Mesh(unit_cube.GenerateMesh(maxh=0.6)),
    "hex": lambda: MakeStructured3DMesh(hexes=True, nx=2, ny=2, nz=2),
    "prism": lambda: MakeStructured3DMesh(hexes=False, prism=True, nx=2, ny=2, nz=2),
}

def tang(v, n):
    return v - (v*n)*n

spaces = {
    "HDiv": (lambda mesh, k, **kw: HDiv(mesh, order=k, **kw),
             lambda s, mu, n: (s*n)*mu, "FacetFESpace"),
    "HDivDiv": (lambda mesh, k, **kw: HDivDiv(mesh, order=k, **kw),
                lambda s, mu, n: ((s*n)*n)*(mu*n), "NormalFacetFESpace"),
    "HCurlDiv": (lambda mesh, k, **kw: HCurlDiv(mesh, order=k, **kw),
                 lambda s, mu, n: tang(s*n, n)*tang(mu, n), "TangentialFacetFESpace"),
}

# multiplier is exact: trace map onto it is surjective, and its kernel is the conforming space
@pytest.mark.parametrize("meshname", meshes)
@pytest.mark.parametrize("spacename", spaces)
@pytest.mark.parametrize("order", [0, 1, 2])
def test_multiplier_exact(meshname, spacename, order):
    if meshname == "prism" and spacename == "HCurlDiv":
        pytest.skip("HCurlDiv has no prism element")
    mesh = meshes[meshname]()
    create, trace, mname = spaces[spacename]
    Xdc = Discontinuous(create(mesh, order))
    M = Xdc.CreateMultiplierSpace()
    assert type(M).__name__ == mname
    s, mu = Xdc.TrialFunction(), M.TestFunction()
    n = specialcf.normal(mesh.dim)
    b = BilinearForm(trialspace=Xdc, testspace=M)
    b += trace(s, mu, n)*dx(element_boundary=True)
    b.Assemble()
    sv = np.linalg.svd(b.mat.ToDense().NumPy(), compute_uv=False)
    rank = np.sum(sv > 1e-10*sv[0])
    assert rank == M.ndof
    Xc = create(mesh, order, dirichlet=".*")
    assert Xdc.ndof - rank == sum(Xc.FreeDofs())

@pytest.mark.parametrize("quad", [False, True])
@pytest.mark.parametrize("order", [0, 1, 2])
def test_multiplier_exact_surface(quad, order):
    from netgen.occ import Box, OCCGeometry
    from netgen.meshing import MeshingStep
    mesh = Mesh(OCCGeometry(Box((0,0,0), (1,1,1))).GenerateMesh(maxh=0.5, quad_dominated=quad,
                                                               perfstepsend=MeshingStep.MESHSURFACE))
    Xdc = Discontinuous(HDivSurface(mesh, order=order), BND=True)
    M = Xdc.CreateMultiplierSpace()
    assert type(M).__name__ == "FacetSurface"
    s, mu = Xdc.TrialFunction(), M.TestFunction()
    nu = Cross(specialcf.tangential(3), specialcf.normal(3))
    b = BilinearForm(trialspace=Xdc, testspace=M)
    b += (s.Trace()*nu)*mu*ds(element_boundary=True)
    b.Assemble()
    sv = np.linalg.svd(b.mat.ToDense().NumPy(), compute_uv=False)
    rank = np.sum(sv > 1e-10*sv[0])
    assert rank == M.ndof
    assert Xdc.ndof - rank == HDivSurface(mesh, order=order).ndof

# SurfaceL2 from HDivSurface.CreateDivergenceSpace contains div and is onto
@pytest.mark.parametrize("quad", [False, True])
@pytest.mark.parametrize("RT", [False, True])
@pytest.mark.parametrize("hodivfree", [False, True])
@pytest.mark.parametrize("order", [0, 1, 2, 3])
def test_surface_divergence_space(quad, RT, hodivfree, order):
    from netgen.occ import Box, OCCGeometry
    from netgen.meshing import MeshingStep
    mesh = Mesh(OCCGeometry(Box((0,0,0), (1,1,1))).GenerateMesh(maxh=0.5, quad_dominated=quad,
                                                               perfstepsend=MeshingStep.MESHSURFACE))
    fes = HDivSurface(mesh, order=order, RT=RT, hodivfree=hodivfree)
    Q = fes.CreateDivergenceSpace()
    X = Discontinuous(fes, BND=True)
    u, v = X.TnT()
    p, q = Q.TnT()
    D = BilinearForm(div(u.Trace())*div(v.Trace())*ds).Assemble().mat.ToDense().NumPy()
    M = BilinearForm(p*q*ds).Assemble().mat.ToDense().NumPy()
    B = BilinearForm(trialspace=X, testspace=Q)
    B += div(u.Trace())*q*ds
    B = B.Assemble().mat.ToDense().NumPy()
    assert np.linalg.norm(D - B.T @ np.linalg.solve(M, B)) < 1e-8*np.linalg.norm(D)
    sv = np.linalg.svd(B, compute_uv=False)
    assert np.sum(sv > 1e-10*sv[0]) == Q.ndof

def test_multiplier_only_on_discontinuous():
    mesh = meshes["trig"]()
    assert not hasattr(HDiv(mesh, order=2), "CreateMultiplierSpace")

def test_multiplier_flags():
    mesh = meshes["trig"]()
    M = Discontinuous(HDiv(mesh, order=2, complex=True)).CreateMultiplierSpace(dirichlet="left|right")
    assert M.is_complex
    assert sum(M.FreeDofs()) < M.ndof
    assert M.globalorder == 2
    with pytest.raises(Exception):
        Discontinuous(H1(mesh, order=2)).CreateMultiplierSpace()

# hybrid mixed Poisson reproduces the conforming mixed solution
def test_hybrid_mixed_poisson():
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.2))
    V = Discontinuous(HDiv(mesh, order=2))
    Q = L2(mesh, order=1)
    M = V.CreateMultiplierSpace(dirichlet=".*")
    X = V*Q*M
    (s, p, l), (ts, q, tl) = X.TnT()
    n = specialcf.normal(2)
    a = BilinearForm((s*ts + div(s)*q + div(ts)*p)*dx
                     - ((s*n)*tl + (ts*n)*l)*dx(element_boundary=True)).Assemble()
    f = LinearForm(-q*dx).Assemble()
    gf = GridFunction(X)
    gf.vec.data = a.mat.Inverse(X.FreeDofs(), inverse="umfpack") * f.vec

    X2 = HDiv(mesh, order=2)*Q
    (s, p), (ts, q) = X2.TnT()
    a2 = BilinearForm((s*ts + div(s)*q + div(ts)*p)*dx).Assemble()
    f2 = LinearForm(-q*dx).Assemble()
    gf2 = GridFunction(X2)
    gf2.vec.data = a2.mat.Inverse(inverse="umfpack") * f2.vec
    assert sqrt(Integrate((gf.components[1]-gf2.components[1])**2, mesh)) < 1e-10

# NormalFacet boundary data: Interpolate == Set (facet-wise L2 projection of u.n), exact for polynomial u.n
@pytest.mark.parametrize("meshname", meshes)
@pytest.mark.parametrize("order", [0, 1, 2])
def test_normalfacet_boundary_interpolate(meshname, order):
    mesh = meshes[meshname]()
    d = mesh.dim
    n = specialcf.normal(d)
    bnd = mesh.Boundaries(".*")
    cf = CF((x**order+2*y, x-y**order) if d == 2 else (x**order+2*y, x-y**order, z**order*x+1))
    H = NormalFacetFESpace(mesh, order=order)
    a, b = GridFunction(H), GridFunction(H)
    a.Set(cf, definedon=bnd)
    b.Interpolate(cf, definedon=bnd)
    assert Norm(a.vec - b.vec) < 1e-12
    if order > 0:
        assert sqrt(Integrate((b*n - cf*n)**2, mesh, BND, order=2*order+4)) < 1e-12

# dual-moment interpolation reproduces boundary FE functions, also on non-affine (mapped) facets
@pytest.mark.parametrize("meshname", ["quad", "hex", "prism"])
@pytest.mark.parametrize("order", [1, 3])
def test_normalfacet_interpolate_reproduces(meshname, order):
    warp = {"quad": lambda x, y: (x+0.1*x*y, y+0.05*x*x),
            "hex": lambda x, y, z: (x+0.1*y*z, y+0.05*x*x, z+0.1*x*y),
            "prism": lambda x, y, z: (x+0.1*y*z, y, z+0.1*x*y)}[meshname]
    mesh = (MakeStructured2DMesh(quads=True, nx=3, ny=3, mapping=warp) if meshname == "quad" else
            MakeStructured3DMesh(hexes=(meshname == "hex"), prism=(meshname == "prism"), nx=2, ny=2, nz=2, mapping=warp))
    bnd = mesh.Boundaries(".*")
    H = NormalFacetFESpace(mesh, order=order)
    gf, gi = GridFunction(H), GridFunction(H)
    mask = H.GetDofs(bnd)
    gf.vec.FV().NumPy()[:] = np.random.RandomState(1).uniform(-1, 1, H.ndof) * np.array(mask)
    gi.Interpolate(gf, definedon=bnd)
    assert np.abs((gi.vec.FV().NumPy() - gf.vec.FV().NumPy())[np.array(mask)]).max() < 1e-12
