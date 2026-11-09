import pytest

from netgen.csg import *
from netgen.geom2d import SplineGeometry
from ngsolve import *
from numpy import pi


def test_periodicH1():
    # incoming wave direction
    d=[1./sqrt(3),1./sqrt(3),1/sqrt(3)]
    
    # wave number
    k = 10
    Lx = 2*pi/(k*d[0])
    Ly = 2*pi/(k*d[1])
    Lz = 2*pi/(k*d[2])
    print(Lx)

    geo = CSGeometry()

    left = Plane(Pnt(0,0,0),Vec(-1,0,0))
    right = Plane(Pnt(Lx,0,0),Vec(1,0,0))
    back = Plane(Pnt(0,0,0),Vec(0,-1,0))
    front = Plane(Pnt(0,Ly,0),Vec(0,1,0))
    bot = Plane(Pnt(0,0,0),Vec(0,0,-1))
    top = Plane(Pnt(0,0,Lz),Vec(0,0,1))

    
    cube = left * right * top * bot * back * front
    
    geo.Add(cube)
    geo.PeriodicSurfaces(left,right)
    geo.PeriodicSurfaces(back,front)
    geo.PeriodicSurfaces(bot,top)

    mesh = Mesh(geo.GenerateMesh(maxh=min(Lx,Ly,Lz)/10))
    
    fes = Periodic(H1(mesh,order=2, complex=True))
    
    u,v = fes.TrialFunction(), fes.TestFunction()
    
    a = BilinearForm(fes)
    a += SymbolicBFI(grad(u)*grad(v) + u*v)

    c = Preconditioner(a,"direct")

    f = LinearForm(fes)
    f += SymbolicLFI((k*k+1) * exp(1J * k * (d[0] * x + d[1] * y + d[2] * z)) * v)

    u = GridFunction(fes)

    u_exact = exp(1J*k*(d[0] * x + d[1] * y + d[2] * z))

    with TaskManager():
        a.Assemble()
        f.Assemble()
        ainv = CGSolver(a.mat,c.mat)
        u.vec.data = ainv * f.vec
        error = sqrt(Integrate(Conj(u-u_exact)*(u-u_exact),mesh).real)
        assert error < 1e-2


def test_periodic_set_order():
    geometry = SplineGeometry()
    points = [geometry.AppendPoint(*point)
              for point in ((0, 0), (1, 0), (1, 1), (0, 1))]
    geometry.Append(["line", points[0], points[1]], bc="outer")
    master = geometry.Append(["line", points[1], points[2]], bc="periodic")
    geometry.Append(["line", points[2], points[3]], bc="outer")
    geometry.Append(["line", points[0], points[3]], leftdomain=0,
                    rightdomain=1, copy=master, bc="periodic")

    mesh = Mesh(geometry.GenerateMesh(maxh=0.2))
    base_space = TangentialFacetFESpace(
        mesh, order=0, order_policy=ORDER_POLICY.VARIABLE)
    fes = Periodic(base_space)

    periodic_edges = mesh.GetPeriodicNodePairs(EDGE)
    assert periodic_edges
    (master_edge, minion_edge), _ = periodic_edges[0]

    fes.SetOrder(NodeId(EDGE, master_edge), 1)
    assert base_space.GetOrder(NodeId(EDGE, master_edge)) == 1
    assert base_space.GetOrder(NodeId(EDGE, minion_edge)) == 1

    fes.SetOrder(NodeId(EDGE, minion_edge), 2)
    assert base_space.GetOrder(NodeId(EDGE, master_edge)) == 2
    assert base_space.GetOrder(NodeId(EDGE, minion_edge)) == 2


def test_quasiperiodic():
    k=1.
    Lx=1.
    Ly=2.
    Lz=1.5
    d=[1./sqrt(3),1./sqrt(3),1./sqrt(3)]

    geo = CSGeometry()

    left = Plane(Pnt(0,0,0),Vec(-1,0,0))
    right = Plane(Pnt(Lx,0,0),Vec(1,0,0))
    bot = Plane(Pnt(0,0,0),Vec(0,0,-1))
    top = Plane(Pnt(0,0,Lz),Vec(0,0,1))
    back = Plane(Pnt(0,0,0),Vec(0,-1,0))
    front = Plane(Pnt(0,Ly,0),Vec(0,1,0))

    cube = left * right * top * bot * back * front

    geo.Add(cube)
    geo.PeriodicSurfaces(left,right)
    geo.PeriodicSurfaces(back,front)
    geo.PeriodicSurfaces(bot,top)

    mesh = Mesh(geo.GenerateMesh(maxh=0.3))
    fx = exp(1J*k*d[0]*Lx)
    fy = exp(1J*k*d[1]*Ly)
    fz = exp(1J*k*d[2]*Lz)
    factors = [fx,fy,fz]
    fes = Periodic(H1(mesh,order=5,complex=True),phase=factors)
    u,v = fes.TnT()

    a = BilinearForm(fes)
    a += SymbolicBFI(grad(u) * grad(v) + u * v)
    f = LinearForm(fes)
    f += SymbolicLFI((k**2+1)*exp(1J*k*(d[0]*x+d[1]*y+d[2]*z)) * v)

    u = GridFunction(fes,name="u")

    with TaskManager():
        a.Assemble()
        f.Assemble()
        u.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec

    u_true = exp(1J * k * (d[0] * x + d[1] * y + d[2] * z))
    assert sqrt(Integrate(Conj(u_true-u)*(u_true-u),mesh).real) < 1e-8


_dual_set_meshes = {}

def _dual_set_mesh(dim):
    """unit-ish square / box, every pair of opposite sides identified periodically"""
    if dim not in _dual_set_meshes:
        import netgen.occ as occ
        per = occ.IdentificationType.PERIODIC
        if dim == 2:
            shape = occ.Rectangle(1, 1.2).Face()
            sides, maxh = shape.edges, 0.3
            dirs = (occ.X, occ.Y)
        else:
            shape = occ.Box((0, 0, 0), (1, 1.2, 1.5))
            sides, maxh = shape.faces, 0.4
            dirs = (occ.X, occ.Y, occ.Z)
        for d, name in zip(dirs, ("px", "py", "pz")):
            sides.Min(d).Identify(sides.Max(d), name, per)
        mesh = Mesh(occ.OCCGeometry(shape, dim=dim).GenerateMesh(maxh=maxh))
        # no element may contain a vertex together with its periodic image
        fes = Periodic(H1(mesh, order=1))
        for el in fes.Elements(VOL):
            assert len(set(el.dofs)) == len(el.dofs)
        _dual_set_meshes[dim] = mesh
    return _dual_set_meshes[dim]


@pytest.mark.parametrize("dim", [2, 3])
@pytest.mark.parametrize("space", [H1, HCurl, HDiv, HCurlCurl, HDivDiv])
def test_periodic_dual_set(space, dim):
    # The dual interpolant Set(cf, dual=True) is element-local, and for a periodic cf the
    # identified dofs get the same value from both sides: Periodic(V) gives the interpolant of V.
    # Spaces caching the inverse duality matrix per vertex-order class must take the class from
    # the vertex numbers PeriodicFESpace::GetFE sets, not from the mesh numbering.
    mesh = _dual_set_mesh(dim)
    s = [sin(2*pi*x), cos(2*pi*y/1.2), sin(2*pi*z/1.5)][:dim]
    if space is H1:
        cf = s[0] * s[1] * (s[2] if dim == 3 else 1) + 0.3
    elif space in (HCurl, HDiv):
        cf = CF(tuple(s[i] * s[(i+1) % dim] + 0.1*(i+1) for i in range(dim)))
    else:
        m = CF(tuple(s[i] * s[j] + 0.2*(i == j) + 0.1*j for i in range(dim) for j in range(dim)), dims=(dim, dim))
        cf = 0.5 * (m + m.trans)
    for order in range(1 if space is H1 else 0, 4):
        gper = GridFunction(Periodic(space(mesh, order=order)))
        gper.Set(cf, dual=True)
        g = GridFunction(space(mesh, order=order))
        g.Set(cf, dual=True)
        err = sqrt(Integrate(InnerProduct(gper - g, gper - g), mesh))
        assert err < 1e-10, f"{space.__name__} order {order}: |periodic - plain| = {err}"
