"""Order specification of H1: SetOrder for element types, regions, elements and nodes,
the minimum rule, orders with a direction, hp-refined meshes, pickling."""
import pickle
from math import cos, sin, floor
import numpy as np
import pytest
from netgen.occ import unit_square, unit_cube, Rectangle, Box, OCCGeometry, Glue, X, Y, Z
from ngsolve import *
from ngsolve.meshes import MakeStructured2DMesh, MakeStructured3DMesh


# -------------------- orders on regions

def consistent(fes):
    return all(fes.GetFE(el).ndof == len(el.dofs) for vb in [VOL, BND] for el in fes.Elements(vb))


def mesh2d(quads=False, maxh=0.3):
    r1 = Rectangle(1, 1).Face(); r1.faces.name = "a"
    r2 = Rectangle(1, 1).Face().Move((1, 0, 0)); r2.faces.name = "b"
    shape = Glue([r1, r2])
    shape.edges.Min(Y).name = "bot"
    shape.edges.Max(Y).name = "top"
    shape.edges.Min(X).name = "left"
    return Mesh(OCCGeometry(shape, dim=2).GenerateMesh(maxh=maxh, quad_dominated=quads))


def mesh3d(maxh=0.6):
    b1 = Box((0, 0, 0), (1, 1, 1)); b1.mat("a")
    b2 = Box((1, 0, 0), (2, 1, 1)); b2.mat("b")
    shape = Glue([b1, b2])
    shape.faces.Min(Z).name = "bot"
    shape.faces.Min(Y).name = "front"
    shape.edges.Min(Y + Z).name = "rim"
    return Mesh(OCCGeometry(shape).GenerateMesh(maxh=maxh))


def ninner(et, p):
    return {ET.SEGM: max(p-1, 0), ET.TRIG: (p-1)*(p-2)//2 if p > 2 else 0, ET.QUAD: (p-1)**2 if p > 1 else 0,
            ET.TET: (p-1)*(p-2)*(p-3)//6 if p > 3 else 0}[et]


def expected(mesh, pcell, bnd_edges={}, bnd_faces={}):
    """orders by the minimum rule, and the number of dofs"""
    pe = {}; pf = {}
    for el in mesh.Elements(VOL):
        p = pcell(el)
        for e in el.edges: pe[e.nr] = min(pe.get(e.nr, 99), p)
        if mesh.dim == 3:
            for f in el.faces: pf[f.nr] = min(pf.get(f.nr, 99), p)
    pe.update(bnd_edges); pf.update(bnd_faces)
    nd = mesh.nv + sum(p-1 for p in pe.values())
    nd += sum(ninner(ET.TRIG, p) for p in pf.values())
    nd += sum(ninner(el.type, pcell(el)) for el in mesh.Elements(VOL))
    return pe, pf, nd


def jumps(fes):
    """squared jumps of a function with arbitrary coefficients over the interior facets"""
    u, v = fes.TnT()
    a = BilinearForm(fes, nonassemble=True)
    a += (u-u.Other())*(v-v.Other())*dx(skeleton=True)
    gf = GridFunction(fes)
    gf.vec.FV().NumPy()[:] = [((7*i) % 11)/11 - 0.4 for i in range(fes.ndof)]
    res = gf.vec.CreateVector()
    a.Apply(gf.vec, res)
    return abs(InnerProduct(res, gf.vec))


def conforming(fes, mesh):
    return jumps(fes) < 1e-20


def test_jump_control():
    mesh = mesh2d()
    assert jumps(H1(mesh, order=3)) < 1e-20
    assert jumps(Discontinuous(H1(mesh, order=3))) > 1e-3
    assert jumps(L2(mesh, order=3)) > 1e-3


@pytest.mark.parametrize("quads", [False, True])
def test_volume_region_2d(quads):
    mesh = mesh2d(quads)
    fes = H1(mesh, order=2)
    fes.SetOrder(mesh.Materials("b"), 4)
    fes.Update()
    pe, pf, nd = expected(mesh, lambda el: 4 if el.mat == "b" else 2)
    assert fes.ndof == nd
    assert all(fes.GetOrder(NodeId(EDGE, i)) == p for i, p in pe.items())
    assert set(pe.values()) == {2, 4}
    assert consistent(fes) and conforming(fes, mesh)
    assert fes.order_policy == ORDER_POLICY.VARIABLE
    # contains the global order
    gf = GridFunction(fes)
    gf.Set(x*x - 3*x*y)
    assert Integrate((gf - (x*x - 3*x*y))**2 * dx, mesh) < 1e-24

    # kept by refinement
    mesh.Refine()
    fes.Update()
    assert fes.ndof == expected(mesh, lambda el: 4 if el.mat == "b" else 2)[2]
    assert consistent(fes) and conforming(fes, mesh)


def test_volume_region_3d():
    mesh = mesh3d()
    fes = H1(mesh, order=2)
    fes.SetOrder(mesh.Materials("a"), 5)
    fes.Update()
    pe, pf, nd = expected(mesh, lambda el: 5 if el.mat == "a" else 2)
    assert fes.ndof == nd
    assert all(fes.GetOrder(NodeId(EDGE, i)) == p for i, p in pe.items())
    assert all(fes.GetOrder(NodeId(FACE, i)) == p for i, p in pf.items())
    assert consistent(fes) and conforming(fes, mesh)


def test_later_call_wins():
    mesh = mesh2d()
    fes = H1(mesh, order=2)
    fes.SetOrder(mesh.Materials("b"), 4)
    fes.SetOrder(mesh.Materials(".*"), 3)
    fes.Update()
    assert fes.ndof == H1(mesh, order=3).ndof
    fes.SetOrder(mesh.Materials("b"), 5)
    fes.SetOrder(mesh.Materials("b"), 4)
    fes.Update()
    assert fes.ndof == expected(mesh, lambda el: 4 if el.mat == "b" else 3)[2]


def test_boundary_region_2d():
    mesh = mesh2d()
    fes = H1(mesh, order=2)
    fes.SetOrder(mesh.Boundaries("bot"), 5)
    fes.Update()
    bot = [el.edges[0].nr for el in mesh.Elements(BND) if el.mat == "bot"]
    assert fes.ndof == H1(mesh, order=2).ndof + 3*len(bot)
    assert all(fes.GetOrder(NodeId(EDGE, i)) == (5 if i in bot else 2) for i in range(mesh.nedge))
    assert consistent(fes) and conforming(fes, mesh)
    # the trace on the boundary has the high order
    gf = GridFunction(fes)
    gf.Set(x**5, definedon=mesh.Boundaries("bot"))
    assert Integrate((gf - x**5)**2 * ds("bot"), mesh) < 1e-24
    # two boundary regions: the later one gets the common vertex edges, here they share none; lowering works too
    fes.SetOrder(mesh.Boundaries("bot|top"), 1)
    fes.Update()
    top = [el.edges[0].nr for el in mesh.Elements(BND) if el.mat == "top"]
    assert fes.ndof == H1(mesh, order=2).ndof - len(bot) - len(top)
    assert consistent(fes) and conforming(fes, mesh)


def test_boundary_region_3d():
    mesh = mesh3d()
    base = H1(mesh, order=2).ndof
    fes = H1(mesh, order=2)
    fes.SetOrder(mesh.Boundaries("bot"), 4)
    fes.Update()
    faces = {mesh[el].faces[0].nr for el in mesh.Elements(BND) if el.mat == "bot"}
    edges = {e.nr for el in mesh.Elements(BND) if el.mat == "bot" for e in el.edges}
    assert fes.ndof == base + 2*len(edges) + 3*len(faces)
    assert all(fes.GetOrder(NodeId(EDGE, i)) == (4 if i in edges else 2) for i in range(mesh.nedge))
    assert all(fes.GetOrder(NodeId(FACE, i)) == (4 if i in faces else 2) for i in range(mesh.nface))
    assert consistent(fes) and conforming(fes, mesh)

    # an edge region, and the last region wins on common edges
    rim = {e.nr for el in mesh.Elements(BBND) if el.mat == "rim" for e in el.edges}
    assert rim and rim <= edges
    fes.SetOrder(mesh.BBoundaries("rim"), 6)
    fes.Update()
    assert fes.ndof == base + 2*len(edges) + 3*len(faces) + 2*len(rim)
    fes.SetOrder(mesh.Boundaries("bot"), 4)        # again: now later than the rim
    fes.Update()
    assert fes.ndof == base + 2*len(edges) + 3*len(faces)
    front_edges = {e.nr for el in mesh.Elements(BND) if el.mat == "front" for e in el.edges}
    fes.SetOrder(mesh.Boundaries("front"), 3)
    fes.Update()
    assert all(fes.GetOrder(NodeId(EDGE, i)) == (3 if i in front_edges else 4 if i in edges else 2) for i in range(mesh.nedge))
    assert consistent(fes) and conforming(fes, mesh)


def test_with_other_layers():
    mesh = mesh2d()
    fes = H1(mesh, order=2, orderinner=3)
    fes.SetOrder(mesh.Materials("b"), 4)
    fes.SetOrder(TRIG, 3)          # layer 5, overridden by orderinner
    fes.SetOrder(NodeId(EDGE, 0), 6)
    fes.Update()
    pe, pf, nd = expected(mesh, lambda el: 4 if el.mat == "b" else 2)
    nd += sum(ninner(ET.TRIG, 3) - ninner(ET.TRIG, 4 if el.mat == "b" else 2) for el in mesh.Elements(VOL))
    nd += 6 - pe[0]
    assert fes.ndof == nd and consistent(fes) and conforming(fes, mesh)

    fes = H1(mesh, order=2, trigorder=3)
    fes.SetOrder(mesh.Materials("b"), 4)
    fes.Update()
    assert fes.ndof == sum(ninner(ET.TRIG, 5 if el.mat == "b" else 3) for el in mesh.Elements(VOL)) \
        + mesh.nv + sum(p-1 for p in expected(mesh, lambda el: 4 if el.mat == "b" else 2)[0].values())
    assert consistent(fes)


def test_vectorh1_and_errors():
    mesh = mesh2d()
    V = VectorH1(mesh, order=2)
    V.SetOrder(mesh.Materials("b"), 4)
    with pytest.raises(Exception, match="order was changed"):
        V.GetFE(ElementId(VOL, 0))
    V.Update()
    assert V.ndof == 2*expected(mesh, lambda el: 4 if el.mat == "b" else 2)[2]
    assert all(V.GetFE(el).ndof == len(el.dofs) for el in V.Elements(VOL))

    with pytest.raises(Exception, match="not available"):
        L2(mesh, order=2).SetOrder(mesh.Materials("b"), 4)
    with pytest.raises(Exception, match="different mesh"):
        H1(mesh, order=2).SetOrder(mesh2d().Materials("b"), 4)
    with pytest.raises(Exception, match="variableorder"):
        H1(mesh, variableorder=True).SetOrder(mesh.Materials("b"), 4)
    fes = H1(mesh, variableorder=True, relorder=1)
    fes.SetOrder(mesh.Boundaries("bot"), 4)
    fes.Update()
    assert consistent(fes)


def test_definedon():
    mesh = mesh2d()
    fes = H1(mesh, order=3, definedon="b")
    fes.SetOrder(mesh.Materials("a"), 1)       # not part of the space: no effect on the interface
    fes.Update()
    assert sum(fes.FreeDofs()) == sum(H1(mesh, order=3, definedon="b").FreeDofs())
    assert consistent(fes)


# -------------------- sequence of SetOrder calls, order policy, stale spaces

def orders(fes, mesh):
    return ([fes.GetOrder(NodeId(EDGE, i)) for i in range(mesh.nedge)],
            [fes.GetOrder(NodeId(FACE, i)) for i in range(mesh.nface)],
            [tuple(el.dofs) for el in fes.Elements(VOL)])


def meshes():
    return [Mesh(unit_square.GenerateMesh(maxh=0.4)), MakeStructured2DMesh(quads=True, nx=2, ny=2),
            Mesh(unit_cube.GenerateMesh(maxh=0.7))]


@pytest.mark.parametrize("space", [H1, VectorH1])
def test_setorder_any_sequence(space):
    for mesh in meshes():
        et = TRIG if mesh.dim == 2 else TET
        def nodes(fes):
            for e in mesh.edges:
                if e.nr % 2: fes.SetOrder(e, 4)
        a = space(mesh, order=2); a.SetOrder(et, 3); nodes(a); a.Update()
        b = space(mesh, order=2); nodes(b); b.SetOrder(et, 3); b.Update()
        c = space(mesh, order=2); nodes(c); c.Update(); c.SetOrder(et, 3); c.Update()
        assert a.ndof == b.ndof == c.ndof
        for f in (a, b, c):
            assert consistent(f)
        if space is H1:
            assert orders(a, mesh) == orders(b, mesh) == orders(c, mesh)
            assert all(a.GetOrder(NodeId(EDGE, i)) == 4 for i in range(1, mesh.nedge, 2))


def test_node_orders_dropped_by_refinement():
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.4))
    fes = H1(mesh, order=2)
    for e in mesh.edges: fes.SetOrder(e, 4)
    fes.Update()
    assert fes.ndof > H1(mesh, order=2).ndof
    assert fes.order_policy == ORDER_POLICY.VARIABLE
    mesh.Refine()
    fes.Update()
    assert fes.ndof == H1(mesh, order=2).ndof
    assert fes.order_policy == ORDER_POLICY.CONSTANT
    # the element type order is kept
    fes.SetOrder(TRIG, 4); fes.Update()
    n = fes.ndof
    mesh.Refine(); fes.Update()
    assert fes.ndof == H1(mesh, order=2, trigorder=4).ndof > n


def test_order_policy_is_derived():
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.4))
    assert H1(mesh, order=3).order_policy == ORDER_POLICY.CONSTANT
    assert H1(mesh, order=3, order_policy=ORDER_POLICY.VARIABLE).order_policy == ORDER_POLICY.CONSTANT
    assert H1(mesh, order=3, orderinner=4).order_policy == ORDER_POLICY.NODETYPE
    assert H1(mesh, order=3, trigorder=4).order_policy == ORDER_POLICY.NODETYPE
    assert H1(mesh, order=2, order_policy=ORDER_POLICY.CONSTANT, orderinner=4).order_policy == ORDER_POLICY.NODETYPE
    fes = H1(mesh, order=2)
    fes.SetOrder(TRIG, 3); fes.Update()
    assert fes.order_policy == ORDER_POLICY.NODETYPE
    fes.SetOrder(NodeId(EDGE, 0), 3); fes.Update()
    assert fes.order_policy == ORDER_POLICY.VARIABLE


def test_constant_policy_flag_is_consistent():
    for mesh in meshes():
        et = TRIG if mesh.dim == 2 else TET
        for p in [1, 2, 3, 4]:
            fes = H1(mesh, order=p, order_policy=ORDER_POLICY.CONSTANT, orderinner=p+2)
            assert consistent(fes) and fes.ndof == H1(mesh, order=p, orderinner=p+2).ndof
            fes = H1(mesh, order=p, order_policy=ORDER_POLICY.CONSTANT)
            fes.SetOrder(et, p+1); fes.Update()
            assert consistent(fes)
            fes = H1(mesh, order=p, order_policy=ORDER_POLICY.CONSTANT)
            fes.SetOrder(NodeId(EDGE, 0), p+1); fes.Update()
            assert consistent(fes) and fes.ndof == H1(mesh, order=p).ndof + 1


def test_setorder_without_update_is_an_error():
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.4))
    for space in [H1, VectorH1, L2]:
        fes = space(mesh, order=2)
        fes.SetOrder(TRIG, 3)
        with pytest.raises(Exception, match="order was changed"):
            fes.GetFE(ElementId(VOL, 0))
        with pytest.raises(Exception, match="order was changed"):
            u, v = fes.TnT()
            BilinearForm(InnerProduct(u, v)*dx).Assemble()
        fes.Update()
        fes.GetFE(ElementId(VOL, 0))
        u, v = fes.TnT()
        BilinearForm(InnerProduct(u, v)*dx).Assemble()

    fes = H1(mesh, order=2)
    fes.SetOrder(NodeId(EDGE, 0), 3)
    with pytest.raises(Exception, match="order was changed"):
        list(fes.Elements(VOL))
    fes.UpdateDofTables()          # enough for single nodes
    assert consistent(fes) and fes.ndof == H1(mesh, order=2).ndof + 1
    fes.SetOrder(TRIG, 3)
    fes.UpdateDofTables()          # not enough for element types
    with pytest.raises(Exception, match="order was changed"):
        list(fes.Elements(VOL))
    fes.Update()
    assert consistent(fes) and fes.ndof == H1(mesh, order=2, trigorder=3).ndof + 1


def test_default_space_solves_like_before():
    # default order policy now takes the fixed order elements
    mesh = Mesh(unit_cube.GenerateMesh(maxh=0.4))
    res = []
    for kw in [dict(), dict(orderinner=3)]:
        fes = H1(mesh, order=3, dirichlet=".*", **kw)
        u, v = fes.TnT()
        a = BilinearForm(grad(u)*grad(v)*dx).Assemble()
        f = LinearForm(x*y*v*dx).Assemble()
        gf = GridFunction(fes)
        gf.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec
        res.append(Integrate(gf*gf*dx, mesh))
    assert abs(res[0]-res[1]) < 1e-14 * abs(res[0])


# -------------------- element orders, inheritance, minimum rule

def pattern(el):
    return 2 + (el.nr * 7) % 4


def set_pattern(fes, mesh):
    for el in mesh.Elements(VOL):
        fes.SetOrder(el, pattern(el))


@pytest.mark.parametrize("make", [lambda: Mesh(unit_square.GenerateMesh(maxh=0.3)),
                                  lambda: Mesh(unit_square.GenerateMesh(maxh=0.3, quad_dominated=True)),
                                  lambda: Mesh(unit_cube.GenerateMesh(maxh=0.5))])
def test_element_orders_minimum_rule(make):
    mesh = make()
    fes = H1(mesh, order=3)
    assert fes.GetOrder(ElementId(VOL, 0)) == 3
    set_pattern(fes, mesh)
    with pytest.raises(Exception, match="order was changed"):
        fes.GetFE(ElementId(VOL, 0))
    fes.Update()
    pe, pf, nd = expected(mesh, pattern)
    assert fes.ndof == nd
    assert all(fes.GetOrder(NodeId(EDGE, i)) == p for i, p in pe.items())
    assert all(fes.GetOrder(NodeId(FACE, i)) == p for i, p in pf.items())
    assert all(fes.GetOrder(el) == pattern(el) for el in mesh.Elements(VOL))
    assert consistent(fes) and conforming(fes, mesh)
    assert fes.order_policy == ORDER_POLICY.VARIABLE
    # an element order is above a region order, a later region does not change it
    fes.SetOrder(mesh.Materials(".*"), 1)
    fes.Update()
    assert fes.ndof == nd


def test_precedence_element_region_global():
    mesh = mesh2d()
    fes = H1(mesh, order=2)
    fes.SetOrder(mesh.Materials("b"), 4)
    el = next(el for el in mesh.Elements(VOL) if el.mat == "b")
    fes.SetOrder(el, 6)
    fes.Update()
    assert [fes.GetOrder(e) for e in mesh.Elements(VOL)] == \
        [6 if e.nr == el.nr else 4 if e.mat == "b" else 2 for e in mesh.Elements(VOL)]
    pcell = lambda e: 6 if e.nr == el.nr else 4 if e.mat == "b" else 2
    assert fes.ndof == expected(mesh, pcell)[2] and consistent(fes) and conforming(fes, mesh)


def inside(V, pt):
    A = np.vstack([np.array(V).T, np.ones(len(V))])
    return np.linalg.solve(A, np.append(pt, 1)).min() > -1e-10


@pytest.mark.parametrize("make", [lambda: Mesh(unit_square.GenerateMesh(maxh=0.3)),
                                  lambda: Mesh(unit_cube.GenerateMesh(maxh=0.5))])
@pytest.mark.parametrize("updates", ["each", "last"])
def test_children_inherit(make, updates):
    mesh = make()
    fes = H1(mesh, order=2)
    set_pattern(fes, mesh)
    fes.Update()
    old = [([np.array(mesh[v].point) for v in el.vertices], pattern(el)) for el in mesh.Elements(VOL)]
    for it in range(2):
        for el in mesh.Elements(VOL):
            mesh.SetRefinementFlag(el, el.nr % 3 == it)
        mesh.Refine()
        if updates == "each": fes.Update()
    fes.Update()
    assert mesh.ne > 3*len(old)
    for el in mesh.Elements(VOL):
        c = sum(np.array(mesh[v].point) for v in el.vertices) / len(el.vertices)
        orders = {p for V, p in old if inside(V, c)}
        assert orders == {fes.GetOrder(el)}
    pcell = lambda el: fes.GetOrder(el)
    assert fes.ndof == expected(mesh, pcell)[2]
    assert consistent(fes) and conforming(fes, mesh)
    # go on: raise the order of some children
    for el in mesh.Elements(VOL):
        if el.nr % 5 == 0: fes.SetOrder(el, fes.GetOrder(el)+1)
    fes.Update()
    assert consistent(fes) and conforming(fes, mesh)


def test_no_parents_drops(capfd):
    mesh = MakeStructured2DMesh(quads=True, nx=3, ny=3)
    fes = H1(mesh, order=2)
    set_pattern(fes, mesh)
    fes.Update()
    assert fes.ndof > H1(mesh, order=2).ndof
    mesh.Refine()
    fes.Update()
    assert "element orders are dropped" in capfd.readouterr().err
    assert fes.ndof == H1(mesh, order=2).ndof
    assert fes.order_policy == ORDER_POLICY.CONSTANT
    assert consistent(fes)


@pytest.mark.parametrize("make,orders", [
    (lambda: MakeStructured2DMesh(quads=True, nx=3, ny=2), [(2, 5), (4, 1), (3, 3)]),
    (lambda: MakeStructured3DMesh(hexes=True, nx=2, ny=2, nz=2), [(2, 5, 3), (4, 1, 2), (3, 3, 3)]),
    (lambda: MakeStructured3DMesh(hexes=False, prism=True, nx=2, ny=2, nz=2), [(2, 2, 5), (4, 4, 1), (3, 3, 3)])])
def test_local_directions(make, orders):
    # the same orders through the mesh (variableorder) and through SetOrder
    mesh = make()
    fes = H1(mesh, order=1)
    for el in mesh.Elements(VOL):
        o = orders[el.nr % len(orders)]
        fes.SetOrder(el, o)
        mesh.SetElementOrders(el, tuple(o) + (o[-1],)*(3-len(o)))
    fes.Update()
    ref = H1(mesh, variableorder=True, relorder=0)
    assert fes.ndof == ref.ndof > H1(mesh, order=1).ndof
    for nt, n in [(EDGE, mesh.nedge), (FACE, mesh.nface)]:
        assert [fes.GetOrder(NodeId(nt, i)) for i in range(n)] == [ref.GetOrder(NodeId(nt, i)) for i in range(n)]
    assert [tuple(el.dofs) for el in fes.Elements(VOL)] == [tuple(el.dofs) for el in ref.Elements(VOL)]
    assert consistent(fes) and conforming(fes, mesh)
    assert consistent(ref) and conforming(ref, mesh)
    # an edge has the lowest order of its cells in its direction
    lo, hi = min(min(o) for o in orders), max(max(o) for o in orders)
    assert all(lo <= fes.GetOrder(NodeId(EDGE, i)) <= hi for i in range(mesh.nedge))
    el0 = ElementId(VOL, 0)
    o0 = orders[0]
    assert fes.GetOrder(el0) == (tuple(o0) + (o0[-1],)*(3-len(o0)) if len(set(o0)) > 1 else o0[0])


def test_mesh_orders_minimum_rule():
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.3))
    for el in mesh.Elements(VOL):
        mesh.SetElementOrder(el, pattern(el))
    fes = H1(mesh, variableorder=True, relorder=1)
    pe, pf, nd = expected(mesh, lambda el: pattern(el)+1)
    assert fes.ndof == nd
    assert all(fes.GetOrder(NodeId(EDGE, i)) == p for i, p in pe.items())
    assert consistent(fes) and conforming(fes, mesh)


def test_hp_refined_mesh():
    face = Rectangle(1, 1).Face()
    face.edges.Min(Y).hpref = 1
    face.vertices.Min(X + Y).hpref = 1
    mesh = Mesh(OCCGeometry(face, dim=2).GenerateMesh(maxh=0.5, quad_dominated=True))
    mesh.RefineHP(levels=3, factor=0.125)
    fes = H1(mesh, variableorder=True, relorder=0)
    assert consistent(fes) and conforming(fes, mesh)
    assert max(fes.GetOrder(NodeId(EDGE, i)) for i in range(mesh.nedge)) == 4


def test_errors_and_vectorh1():
    trig = Mesh(unit_square.GenerateMesh(maxh=0.4))
    fes = H1(trig, order=2)
    with pytest.raises(Exception, match="different orders"):
        fes.SetOrder(ElementId(VOL, 0), (2, 3))
    with pytest.raises(Exception, match="volume elements"):
        fes.SetOrder(ElementId(BND, 0), 3)
    with pytest.raises(Exception, match="not available"):
        L2(trig, order=2).SetOrder(ElementId(VOL, 0), 3)
    prism = MakeStructured3DMesh(hexes=False, prism=True, nx=1, ny=1, nz=1)
    with pytest.raises(Exception, match="prism"):
        H1(prism, order=2).SetOrder(ElementId(VOL, 0), (2, 3, 4))

    V = VectorH1(trig, order=2)
    set_pattern(V, trig)
    V.Update()
    assert V.ndof == 2*expected(trig, pattern)[2]
    assert V.GetOrder(ElementId(VOL, 1)) == pattern(ElementId(VOL, 1))
    assert all(V.GetFE(el).ndof == len(el.dofs) for el in V.Elements(VOL))


# -------------------- orders with a direction, hp-refined meshes

def edge_vec(mesh, e):
    p0, p1 = [np.array(mesh[v].point) for v in e.vertices]
    return p1 - p0, 0.5*(p0+p1)


def along(mesh, e, d):
    t, mid = edge_vec(mesh, e)
    dd = np.array(d(mid) if callable(d) else d, dtype=float)
    return abs(t @ dd) / np.linalg.norm(t) / np.linalg.norm(dd) > 0.7


@pytest.mark.parametrize("direction,pt,pn", [((0, 1), 2, 5), ((1, 0), 2, 5), ((0.2, 1), 4, 1)])
def test_quads_direction(direction, pt, pn):
    mesh = MakeStructured2DMesh(quads=True, nx=3, ny=4)
    fes = H1(mesh, order=3)
    fes.SetOrder(mesh.Materials(".*"), (pt, pn), direction=CF(direction))
    fes.Update()
    pe = [pn if along(mesh, e, direction) else pt for e in mesh.edges]
    assert [fes.GetOrder(NodeId(EDGE, e.nr)) for e in mesh.edges] == pe
    assert fes.ndof == mesh.nv + sum(p-1 for p in pe) + mesh.ne*(pt-1)*(pn-1)
    assert consistent(fes) and conforming(fes, mesh)
    assert all(sorted(set(fes.GetOrder(el))) == sorted({pt, pn}) for el in mesh.Elements(VOL))
    # the same with a tuple as direction, and an order without direction afterwards replaces it
    fes.SetOrder(mesh.Materials(".*"), (pt, pn), direction=direction)
    fes.Update()
    assert [fes.GetOrder(NodeId(EDGE, e.nr)) for e in mesh.edges] == pe
    fes.SetOrder(mesh.Materials(".*"), 3)
    fes.Update()
    assert fes.ndof == H1(mesh, order=3).ndof and fes.GetOrder(ElementId(VOL, 0)) == 3


def test_quads_radial():
    mapping = lambda x, y: ((1+x)*cos(1.2*y), (1+x)*sin(1.2*y))
    mesh = MakeStructured2DMesh(quads=True, nx=3, ny=6, mapping=mapping)
    radial = lambda p: (p[0], p[1])
    for d, tangential in [(CF((x, y)), False), (CF((-y, x)), True)]:
        fes = H1(mesh, order=2)
        fes.SetOrder(mesh.Materials(".*"), (2, 6), direction=d)
        fes.Update()
        pe = [6 if along(mesh, e, radial) != tangential else 2 for e in mesh.edges]
        assert [fes.GetOrder(NodeId(EDGE, e.nr)) for e in mesh.edges] == pe
        assert 2 in pe and 6 in pe
        assert consistent(fes) and conforming(fes, mesh)


def test_hexes_direction():
    mesh = MakeStructured3DMesh(hexes=True, nx=2, ny=3, nz=2)
    for direction in [(0, 0, 1), (1, 0, 0), (0.1, 1, 0.2)]:
        fes = H1(mesh, order=2)
        fes.SetOrder(mesh.Materials(".*"), (2, 5), direction=CF(direction))
        fes.Update()
        axis = int(np.argmax(np.abs(direction)))
        pe = [5 if along(mesh, e, direction) else 2 for e in mesh.edges]
        assert [fes.GetOrder(NodeId(EDGE, e.nr)) for e in mesh.edges] == pe
        # faces across the direction have (2,2), the others (2,5)
        nface = 0
        for f in mesh.faces:
            pts = np.array([mesh[v].point for v in f.vertices])
            nface += 1 if np.ptp(pts[:, axis]) < 1e-10 else 4
        assert fes.ndof == mesh.nv + sum(p-1 for p in pe) + nface + mesh.ne*4
        assert consistent(fes) and conforming(fes, mesh)


def test_prisms_axis():
    mesh = MakeStructured3DMesh(hexes=False, prism=True, nx=2, ny=2, nz=3)
    fes = H1(mesh, order=2)
    fes.SetOrder(mesh.Materials(".*"), (3, 5))
    fes.Update()
    axis_edges = set()
    for el in mesh.Elements(VOL):
        v = [vi.nr for vi in el.vertices]
        for a, b in [(0, 3), (1, 4), (2, 5)]:
            axis_edges |= {e.nr for e in el.edges if {vi.nr for vi in mesh[e].vertices} == {v[a], v[b]}}
    pe = [5 if e.nr in axis_edges else 3 for e in mesh.edges]
    assert [fes.GetOrder(NodeId(EDGE, e.nr)) for e in mesh.edges] == pe
    nface = sum(1 if len(f.vertices) == 3 else 2*4 for f in mesh.faces)
    assert fes.ndof == mesh.nv + sum(p-1 for p in pe) + nface + mesh.ne*4
    assert fes.GetOrder(ElementId(VOL, 0)) == (3, 3, 5)
    assert consistent(fes) and conforming(fes, mesh)


def test_direction_errors():
    trig = Mesh(unit_square.GenerateMesh(maxh=0.4))
    with pytest.raises(Exception, match="only for quadrilaterals"):
        H1(trig, order=2).SetOrder(trig.Materials(".*"), (2, 4), direction=CF((0, 1)))
    quad = MakeStructured2DMesh(quads=True, nx=2, ny=2)
    fes = H1(quad, order=2)
    with pytest.raises(Exception, match="need a direction"):
        fes.SetOrder(quad.Materials(".*"), (2, 4))
    with pytest.raises(Exception, match="one component per space dimension"):
        fes.SetOrder(quad.Materials(".*"), (2, 4), direction=CF((0, 0, 1)))
    with pytest.raises(Exception, match="only for volume regions"):
        fes.SetOrder(quad.Boundaries(".*"), (2, 4), direction=CF((0, 1)))
    with pytest.raises(Exception):
        fes.SetOrder(quad.Materials(".*"), 3, direction=CF((0, 1)))
    with pytest.raises(Exception, match="not available"):
        L2(quad, order=2).SetOrder(quad.Materials(".*"), (2, 4), direction=CF((0, 1)))


# ---------------------------------------------------------------- hp-refined meshes

def hp2d(levels=3):
    face = Rectangle(1, 1).Face()
    face.edges.Min(Y).hpref = 1
    mesh = Mesh(OCCGeometry(face, dim=2).GenerateMesh(maxh=0.5, quad_dominated=True))
    mesh.RefineHP(levels=levels, factor=0.125)
    return mesh


def hp2d_corner(levels=3):
    face = Rectangle(1, 1).Face()
    face.edges.Min(Y).hpref = 1
    face.vertices.Min(X + Y).hpref = 1
    mesh = Mesh(OCCGeometry(face, dim=2).GenerateMesh(maxh=0.5))
    mesh.RefineHP(levels=levels, factor=0.125)
    return mesh


def hp3d(levels=2):
    box = Box((0, 0, 0), (1, 1, 1))
    box.edges.Min(X + Y).hpref = 1
    mesh = Mesh(OCCGeometry(box).GenerateMesh(maxh=0.8))
    mesh.RefineHP(levels=levels, factor=0.125)
    return mesh


def node_orders(fes, mesh):
    return ([fes.GetOrder(NodeId(EDGE, i)) for i in range(mesh.nedge)],
            [fes.GetOrder(NodeId(FACE, i)) for i in range(mesh.nface)])


@pytest.mark.parametrize("make,L", [(hp2d, 3), (hp2d_corner, 3), (hp3d, 2)])
def test_hpslope_is_former_variableorder(make, L):
    mesh = make(L)
    ref = H1(mesh, variableorder=True, relorder=0)
    fes = H1(mesh, order=L+1, hpslope=1)
    assert fes.ndof == ref.ndof > H1(mesh, order=1).ndof
    assert node_orders(fes, mesh) == node_orders(ref, mesh)
    assert [tuple(el.dofs) for el in fes.Elements(VOL)] == [tuple(el.dofs) for el in ref.Elements(VOL)]
    assert consistent(fes) and conforming(fes, mesh)
    assert fes.order_policy == ORDER_POLICY.VARIABLE


@pytest.mark.parametrize("p,slope", [(4, 0.5), (6, 1), (2, 1), (5, 0.34), (4, 2)])
def test_hpslope_formula(p, slope):
    L = 3
    mesh = hp2d(L)
    # the mesh order of an edge is L+1-level
    ref = H1(mesh, variableorder=True, relorder=0)
    fes = H1(mesh, order=p, hpslope=slope)
    used = set(e.nr for el in mesh.Elements(VOL) for e in el.edges)
    for i in used:
        level = L+1 - ref.GetOrder(NodeId(EDGE, i))
        assert fes.GetOrder(NodeId(EDGE, i)) == max(1, p - floor(slope*level + 1e-9))
    assert consistent(fes) and conforming(fes, mesh)
    # along the singular edge the order is not reduced
    par = [i for i in used if abs(edge_vec(mesh, mesh.edges[i])[0][1]) < 1e-12]
    assert par and all(fes.GetOrder(NodeId(EDGE, i)) == p for i in par)


def test_hpslope_with_other_layers():
    mesh = hp2d(3)
    assert H1(mesh, order=3, hpslope=0).ndof == H1(mesh, order=3).ndof
    fes = H1(mesh, order=2, hpslope=1)
    fes.SetOrder(mesh.Materials(".*"), 4)
    fes.Update()
    assert fes.ndof == H1(mesh, order=4, hpslope=1).ndof
    # an element order replaces the rule in its cell
    el = ElementId(VOL, 0)
    fes.SetOrder(el, 4)
    fes.Update()
    assert fes.GetOrder(el) == 4 and consistent(fes) and conforming(fes, mesh)
    with pytest.raises(Exception, match="hpslope"):
        H1(mesh, variableorder=True, hpslope=1)
    V = VectorH1(mesh, order=4, hpslope=1)
    assert V.ndof == 2*H1(mesh, order=4, hpslope=1).ndof
    # no hp-refinement: no levels
    plain = Mesh(unit_square.GenerateMesh(maxh=0.3))
    assert H1(plain, order=3, hpslope=1).ndof == H1(plain, order=3).ndof


# -------------------- pickling

def roundtrip(fes):
    dofs = [tuple(el.dofs) for el in fes.Elements(VOL)]
    f2 = pickle.loads(pickle.dumps(fes))
    assert f2.ndof == fes.ndof and [tuple(el.dofs) for el in f2.Elements(VOL)] == dofs and consistent(f2)
    assert f2.order_policy == fes.order_policy
    f2.Update()
    assert f2.ndof == fes.ndof and [tuple(el.dofs) for el in f2.Elements(VOL)] == dofs
    gf = GridFunction(fes)
    gf.Set(x*y + x*x)
    g2 = pickle.loads(pickle.dumps(gf))
    assert len(g2.vec) == len(gf.vec)
    assert abs(Integrate(g2*g2*dx, g2.space.mesh) - Integrate(gf*gf*dx, fes.mesh)) < 1e-14
    return f2


def cases():
    m = mesh2d()
    q = MakeStructured2DMesh(quads=True, nx=3, ny=3)
    def reg(f): f.SetOrder(m.Materials("b"), 4)
    def bnd(f): f.SetOrder(m.Boundaries("bot"), 4)
    def els(f): [f.SetOrder(el, 2+el.nr % 3) for el in m.Elements(VOL)]
    def nodes(f): [f.SetOrder(e, 4) for e in m.edges]
    def allof(f): reg(f); bnd(f); els(f); nodes(f); f.SetOrder(TRIG, 3)
    def direction(f): f.SetOrder(q.Materials(".*"), (2, 4), direction=CF((0.1, 1)))
    yield "plain", m, lambda f: None
    yield "volume region", m, reg
    yield "boundary region", m, bnd
    yield "element orders", m, els
    yield "node orders", m, nodes
    yield "all layers", m, allof
    yield "direction", q, direction


@pytest.mark.parametrize("space", [H1, VectorH1])
@pytest.mark.parametrize("name", [c[0] for c in cases()])
def test_pickle(space, name):
    mesh, setup = next((c[1], c[2]) for c in cases() if c[0] == name)
    fes = space(mesh, order=2)
    setup(fes)
    fes.Update()
    if space is H1:
        f2 = roundtrip(fes)
    else:
        f2 = pickle.loads(pickle.dumps(fes))
        assert f2.ndof == fes.ndof
        f2.Update()
        assert f2.ndof == fes.ndof
        assert [tuple(el.dofs) for el in f2.Elements(VOL)] == [tuple(el.dofs) for el in fes.Elements(VOL)]
    # a plain space keeps the short state
    assert len(fes.__getstate__()) == (3 if name == "plain" else 4)


def test_unpickled_space_goes_on():
    mesh = mesh2d()
    fes = H1(mesh, order=2)
    fes.SetOrder(mesh.Materials("b"), 4)
    fes.SetOrder(ElementId(VOL, 0), 5)
    fes.Update()
    f2 = pickle.loads(pickle.dumps(fes))
    m2 = f2.mesh
    m2.Refine()
    f2.Update()
    mesh.Refine()
    fes.Update()
    assert f2.ndof == fes.ndof and consistent(f2)
