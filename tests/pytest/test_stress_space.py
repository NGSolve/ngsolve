import pytest
from math import pi, log
from ngsolve import *
from netgen.meshing import Mesh as NM, Element2D, Element1D, MeshPoint, Pnt, FaceDescriptor
from ngsolve.meshes import MakeStructured2DMesh

def square(kind, n):
    # trig, quad, or mixed (every third cell split into two trigs)
    if kind == "quad":
        return MakeStructured2DMesh(quads=True, nx=n, ny=n)
    m = NM(dim=2)
    m.Add(FaceDescriptor(surfnr=1, domin=1, bc=1))
    for i, name in enumerate(["bottom", "right", "top", "left"]):
        m.SetBCName(i, name)
    p = {(i, j): m.Add(MeshPoint(Pnt(i/n, j/n, 0))) for i in range(n+1) for j in range(n+1)}
    for i in range(n):
        for j in range(n):
            a, b, c, d = p[i, j], p[i+1, j], p[i+1, j+1], p[i, j+1]
            if kind == "trig" or (i+j) % 3 == 0:
                m.Add(Element2D(1, [a, b, c])); m.Add(Element2D(1, [a, c, d]))
            else:
                m.Add(Element2D(1, [a, b, c, d]))
    for i in range(n):
        m.Add(Element1D([p[i, 0], p[i+1, 0]], index=1))
        m.Add(Element1D([p[n, i], p[n, i+1]], index=2))
        m.Add(Element1D([p[i+1, n], p[i, n]], index=3))
        m.Add(Element1D([p[0, i+1], p[0, i]], index=4))
    return Mesh(m)

E, nu = 1.0, 0.3
mu, lam = E/2/(1+nu), E*nu/((1+nu)*(1-2*nu))

def tdnns_errors(mesh, k, type1):
    b = sin(pi*x)*sin(pi*y)
    uex = CF((b, 2*x*b))
    gu = CF(tuple(uex[i].Diff(c) for i in range(2) for c in (x, y)), dims=(2, 2))
    sex = mu*(gu+gu.trans) + lam*Trace(gu)*Id(2)
    f = -CF(tuple(sum(sex[i, j].Diff(c) for j, c in enumerate((x, y))) for i in range(2)))
    C = lambda s: 1/(2*mu)*(s - lam/(2*mu+2*lam)*Trace(s)*Id(2))
    Q = HCurl(mesh, order=k, type1=type1, dirichlet=".*")
    X = Q.CreateStressSpace()*Q
    (s, u), (t, v) = X.TnT()
    n = specialcf.normal(2)
    tang = lambda w: w - (w*n)*n
    a = BilinearForm((InnerProduct(C(s), t) + div(s)*v + div(t)*u)*dx(bonus_intorder=4)
                     - ((s*n)*tang(v) + (t*n)*tang(u))*dx(element_boundary=True, bonus_intorder=4)).Assemble()
    l = LinearForm(-f*v*dx).Assemble()
    gf = GridFunction(X)
    gf.vec.data = a.mat.Inverse(X.FreeDofs(), inverse="umfpack")*l.vec
    s_h, u_h = gf.components
    return (sqrt(Integrate(InnerProduct(s_h-sex, s_h-sex), mesh, order=2*k+6)),
            sqrt(Integrate((u_h-uex)**2, mesh, order=2*k+6)))

# sigma converges with rate k, u with rate k (type1) or k+1 (default);
# type1 order 1 (Whitney, piecewise rigid) does not converge with any stress space
@pytest.mark.parametrize("kind", ["trig", "quad", "mixed"])
@pytest.mark.parametrize("type1,k", [(True, 2), (True, 3), (False, 1), (False, 2), (False, 3)])
def test_stress_space_rates(kind, type1, k):
    (es1, eu1), (es2, eu2) = [tdnns_errors(square(kind, n), k, type1) for n in (4, 8)]
    assert log(es1/es2)/log(2) > k - 0.3
    assert log(eu1/eu2)/log(2) > k + (0 if type1 else 1) - 0.4

# stress orders: k-1 (type1, quads), k (default trigs), at least 1 except type1 k=1 without simplices
@pytest.mark.parametrize("kind", ["trig", "quad", "mixed"])
@pytest.mark.parametrize("type1", [True, False])
@pytest.mark.parametrize("k", [1, 2])
def test_stress_space_choice_and_dirichlet(kind, type1, k):
    mesh = square(kind, 3)
    Q = HCurl(mesh, order=k, type1=type1)
    V = Q.CreateStressSpace(dirichlet="top")
    assert type(V).__name__ == "HDivDiv"
    lo = 0 if (type1 and kind == "quad") else 1     # type1 order 1 on quad-only meshes: stress order 0
    for el in mesh.Elements(VOL):
        oi = V.GetOrder(NodeId(ELEMENT, el.nr))
        if type1 or el.type == ET.QUAD:
            assert oi == max(k-1, lo)
        else:
            assert oi == k
    # nn-essential condition: exactly the facet dofs on 'top' are fixed
    top = set(d for el in mesh.Elements(BND) if el.mat == "top" for d in V.GetDofNrs(el))
    fixed = set(i for i, f in enumerate(V.FreeDofs()) if not f)
    assert fixed == top and len(top) > 0
    assert V.ndof < Q.CreateStressSpace(discontinuous=True).ndof

# 3D tets: type1 HDivDiv(max(k-1,1)), default HDivDiv(max(k-1,1), plus)
from ngsolve.meshes import MakeStructured3DMesh

def tdnns_errors_3d(mesh, k, type1):
    b = sin(pi*x)*sin(pi*y)*sin(pi*z)
    uex = CF((b, 2*x*b, (1+y*y)*b))
    gu = CF(tuple(uex[i].Diff(c) for i in range(3) for c in (x, y, z)), dims=(3, 3))
    sex = mu*(gu+gu.trans) + lam*Trace(gu)*Id(3)
    f = -CF(tuple(sum(sex[i, j].Diff(c) for j, c in enumerate((x, y, z))) for i in range(3)))
    C = lambda s: 1/(2*mu)*(s - lam/(2*mu+3*lam)*Trace(s)*Id(3))
    Q = HCurl(mesh, order=k, type1=type1, dirichlet=".*")
    X = Q.CreateStressSpace()*Q
    (s, u), (t, v) = X.TnT()
    n = specialcf.normal(3)
    tang = lambda w: w - (w*n)*n
    a = BilinearForm((InnerProduct(C(s), t) + div(s)*v + div(t)*u)*dx
                     - ((s*n)*tang(v) + (t*n)*tang(u))*dx(element_boundary=True)).Assemble()
    l = LinearForm(-f*v*dx).Assemble()
    gf = GridFunction(X)
    gf.vec.data = a.mat.Inverse(X.FreeDofs(), inverse="umfpack")*l.vec
    s_h, u_h = gf.components
    return (sqrt(Integrate(InnerProduct(s_h-sex, s_h-sex), mesh, order=2*k+4)),
            sqrt(Integrate((u_h-uex)**2, mesh, order=2*k+4)))

# sigma rate k-1 (type1) / k (default), u rate k / k+1
@pytest.mark.parametrize("type1,k", [(True, 2), (False, 1), (False, 2)])
def test_stress_space_tet_rates(type1, k):
    (es1, eu1), (es2, eu2) = [tdnns_errors_3d(MakeStructured3DMesh(hexes=False, nx=n, ny=n, nz=n), k, type1)
                              for n in (2, 4)]
    ks = k-1 if type1 else k
    assert log(es1/es2)/log(2) > ks - 0.5
    assert log(eu1/eu2)/log(2) > ks + 1 - 0.5

@pytest.mark.parametrize("type1", [True, False])
@pytest.mark.parametrize("k", [1, 2, 3])
def test_stress_space_tet_choice(type1, k):
    mesh = MakeStructured3DMesh(hexes=False, nx=2, ny=2, nz=2)
    Q = HCurl(mesh, order=k, type1=type1)
    V = Q.CreateStressSpace(dirichlet="top")
    ref = HDivDiv(mesh, order=max(k-1, 1), plus=not type1, dirichlet="top")
    assert V.ndof == ref.ndof
    assert sum(1 for f in V.FreeDofs() if not f) == sum(1 for f in ref.FreeDofs() if not f) > 0

def pyramid_mesh():
    m = NM(dim=3)
    m.Add(FaceDescriptor(surfnr=1, domin=1, bc=1))
    p = [m.Add(MeshPoint(Pnt(*q))) for q in [(0,0,0), (1,0,0), (1,1,0), (0,1,0), (0.5,0.5,1)]]
    from netgen.meshing import Element3D
    m.Add(Element3D(1, p))
    m.Add(Element2D(1, [p[0], p[3], p[2], p[1]]))
    for a, b in [(0, 1), (1, 2), (2, 3), (3, 0)]:
        m.Add(Element2D(1, [p[a], p[b], p[4]]))
    return Mesh(m)

def test_stress_space_pyramid_raises():
    with pytest.raises(Exception):
        HCurl(pyramid_mesh(), order=2).CreateStressSpace()

def prism_hex_mesh():
    # unit cube cells: one hex, one cell split into two prisms (shared quad faces)
    m = NM(dim=3)
    m.Add(FaceDescriptor(surfnr=1, domin=1, bc=1))
    from netgen.meshing import Element3D
    P = {(i, j, k): m.Add(MeshPoint(Pnt(i, j, k))) for i in range(3) for j in range(2) for k in range(2)}
    m.Add(Element3D(1, [P[0,0,0], P[1,0,0], P[1,1,0], P[0,1,0], P[0,0,1], P[1,0,1], P[1,1,1], P[0,1,1]]))
    m.Add(Element3D(1, [P[1,0,0], P[2,0,0], P[2,1,0], P[1,0,1], P[2,0,1], P[2,1,1]]))
    m.Add(Element3D(1, [P[1,0,0], P[2,1,0], P[1,1,0], P[1,0,1], P[2,1,1], P[1,1,1]]))
    mesh = Mesh(m)
    return mesh

def bl_mesh():
    # tets with a prism boundary layer
    from netgen.occ import Box, OCCGeometry
    from netgen.meshing import BoundaryLayerParameters
    bl = BoundaryLayerParameters(boundary=".*", thickness=[0.1], new_material="layer")
    return Mesh(OCCGeometry(Box((0,0,0), (1,1,1))).GenerateMesh(maxh=0.6, boundary_layers=[bl]))

# nn-conformity on mixed 3D meshes: trace of the broken stress space onto the multiplier is
# surjective and its kernel is the conforming space without the boundary nn-dofs
@pytest.mark.parametrize("meshmaker", [prism_hex_mesh, bl_mesh])
@pytest.mark.parametrize("type1", [True, False])
@pytest.mark.parametrize("k", [1, 2, 3])
def test_stress_space_mixed3d_conforming(meshmaker, type1, k):
    import numpy as np
    mesh = meshmaker()
    V = HCurl(mesh, order=k, type1=type1).CreateStressSpace()
    Xdc = Discontinuous(V)
    M = Xdc.CreateMultiplierSpace()
    n = specialcf.normal(3)
    s, mu = Xdc.TrialFunction(), M.TestFunction()
    b = BilinearForm(trialspace=Xdc, testspace=M)
    b += ((s*n)*n)*(mu*n)*dx(element_boundary=True)
    b.Assemble()
    sv = np.linalg.svd(b.mat.ToDense().NumPy(), compute_uv=False)
    r = int((sv > 1e-10*sv[0]).sum())
    outer = [f for f in range(mesh.nface) if len(mesh[NodeId(FACE, f)].elements) == 1]
    nb = len({d for f in outer for d in V.GetDofNrs(NodeId(FACE, f))})
    assert r == M.ndof
    assert Xdc.ndof - r == V.ndof - nb


# Solve with a scalar Dirichlet value on a vector component of a product space
# (the value used to be interpolated with the wrong dimension); same result as dirichlet flags
def test_solve_component_dirichlet():
    mesh = MakeStructured2DMesh(quads=True, nx=4, ny=2)
    C = lambda s: 1/(2*mu)*(s - lam/(2*mu+2*lam)*Trace(s)*Id(2))
    n = specialcf.normal(2)
    tang = lambda w: w - (w*n)*n
    res = []
    for flags in (False, True):
        Q = HCurl(mesh, order=2, dirichlet="left" if flags else "")
        S = Q.CreateStressSpace(dirichlet="top|bottom|right")
        (s, u), (t, v) = (S*Q).TnT()
        a = (InnerProduct(C(s), t) + div(s)*v + div(t)*u)*dx - ((s*n)*tang(v) + (t*n)*tang(u))*dx(element_boundary=True)
        rhs = -CF((0, -1))*v*dx
        gf = Solve(a == rhs) if flags else Solve(a == rhs, u[BND("left")] == 0)
        res.append(gf.components[1](mesh(1, 0.5)))
    assert abs(res[0][1]-res[1][1]) < 1e-10*abs(res[1][1])
    with pytest.raises(Exception):
        u[BND("left")] == CF((1, 2, 3))

# Solve with automatic static condensation (Discontinuous stress) and Dirichlet data on the
# tangential (HCurl) and normal (NormalFacet multiplier) displacement = the uncondensed solve
@pytest.mark.parametrize("kind", ["trig", "quad", "mixed"])
def test_solve_hybrid_condensed(kind):
    mesh = square(kind, 4)
    C = lambda s: 1/(2*mu)*(s - lam/(2*mu+2*lam)*Trace(s)*Id(2))
    n = specialcf.normal(2)
    tang = lambda w: w - (w*n)*n
    nn = lambda q: (q*n)*n
    uex = CF((sin(x)*y, cos(y)+x))
    f = CF((1, -1))
    Q = HCurl(mesh, order=2)
    S = Discontinuous(Q.CreateStressSpace())
    M = S.CreateMultiplierSpace()
    (s, u, un), (t, v, vn) = (S*Q*M).TnT()
    a = (InnerProduct(C(s), t) + div(s)*v + div(t)*u)*dx \
        - ((s*n)*tang(v) + (t*n)*tang(u) + nn(s)*(vn*n) + nn(t)*(un*n))*dx(element_boundary=True)
    bcs = (u[BND(".*")] == uex, un[BND(".*")] == uex)
    gc = Solve(a == -f*v*dx, *bcs)
    gf = Solve(a == -f*v*dx, *bcs, condense=False)
    assert Norm(gc.vec - gf.vec) < 1e-8 * Norm(gf.vec)
