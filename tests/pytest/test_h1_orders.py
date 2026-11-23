"""Order specification of H1, L2 and HCurl: SetOrder for element types, regions, elements and nodes,
the minimum rule, orders with a direction, hp-refined meshes, pickling, linked spaces."""
import gc
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
        HDiv(mesh, order=2).SetOrder(mesh.Materials("b"), 4)
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
        HDiv(trig, order=2).SetOrder(ElementId(VOL, 0), 3)
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
        HDiv(quad, order=2).SetOrder(quad.Materials(".*"), (2, 4), direction=CF((0, 1)))


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


# -------------------- linked spaces

def edge_orders(fes, mesh):
    used = sorted({e.nr for el in fes.Elements(VOL) for e in el.edges})
    return [fes.GetOrder(NodeId(EDGE, i)) for i in used]


def same_space(a, b):
    return a.ndof == b.ndof and [tuple(el.dofs) for el in a.Elements(VOL)] == [tuple(el.dofs) for el in b.Elements(VOL)]


def test_link_offset():
    mesh = mesh2d()
    V = H1(mesh, order=3)
    W = H1(mesh, orderlink=V, relorder=-1)
    U = H1(mesh, orderlink=V, relorder=2)
    assert same_space(W, H1(mesh, order=2)) and W.globalorder == 2
    assert same_space(U, H1(mesh, order=5)) and U.globalorder == 5
    assert W.order_policy == ORDER_POLICY.CONSTANT

    V.SetOrder(mesh.Materials("b"), 5)
    V.SetOrder(mesh.Boundaries("bot"), 6)
    V.SetOrder(ElementId(VOL, 0), 4)
    V.Update()
    # the companions are not up to date
    for fes in (W, U):
        with pytest.raises(Exception, match="order was changed"):
            fes.GetFE(ElementId(VOL, 0))
        fes.Update()
    ref = H1(mesh, order=2)
    ref.SetOrder(mesh.Materials("b"), 4)
    ref.SetOrder(mesh.Boundaries("bot"), 5)
    ref.SetOrder(ElementId(VOL, 0), 3)
    ref.Update()
    assert same_space(W, ref)
    assert edge_orders(W, mesh) == [p-1 for p in edge_orders(V, mesh)]
    assert edge_orders(U, mesh) == [p+2 for p in edge_orders(V, mesh)]
    assert all(W.GetOrder(el) == V.GetOrder(el)-1 for el in mesh.Elements(VOL))
    for fes in (V, W, U):
        assert consistent(fes) and conforming(fes, mesh)
    assert W.order_policy == ORDER_POLICY.VARIABLE


def test_only_the_creator_sets():
    mesh = mesh2d()
    V = H1(mesh, order=3)
    W = H1(mesh, orderlink=V, relorder=-1)
    with pytest.raises(Exception, match="orders belong to"):
        W.SetOrder(mesh.Materials("b"), 4)
    with pytest.raises(Exception, match="orders belong to"):
        W.SetOrder(ElementId(VOL, 0), 4)
    with pytest.raises(Exception, match="orders belong to"):
        W.SetOrder(mesh.Boundaries("bot"), 4)
    # its own layers
    W.SetOrder(TRIG, 4)
    W.SetOrder(NodeId(EDGE, 0), 5)
    W.Update()
    ref = H1(mesh, order=2, trigorder=4)
    ref.SetOrder(NodeId(EDGE, 0), 5)
    ref.Update()
    assert same_space(W, ref)
    assert same_space(V, H1(mesh, order=3))
    V.SetOrder(mesh.Materials("b"), 5)
    V.Update(); W.Update()
    assert W.GetOrder(NodeId(EDGE, 0)) == 5 and consistent(W) and conforming(W, mesh)


@pytest.mark.parametrize("first", ["creator", "companion"])
def test_refinement_update_order(first):
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.3))
    V = H1(mesh, order=2)
    W = H1(mesh, orderlink=V, relorder=1)
    for el in mesh.Elements(VOL):
        V.SetOrder(el, 2 + el.nr % 3)
    V.Update(); W.Update()
    for el in mesh.Elements(VOL):
        mesh.SetRefinementFlag(el, el.nr % 3 == 0)
    mesh.Refine()
    for fes in ((V, W) if first == "creator" else (W, V)):
        fes.Update()
    assert all(W.GetOrder(el) == V.GetOrder(el)+1 for el in mesh.Elements(VOL))
    assert edge_orders(W, mesh) == [p+1 for p in edge_orders(V, mesh)]
    assert len({V.GetOrder(el) for el in mesh.Elements(VOL)}) == 3
    for fes in (V, W):
        assert consistent(fes) and conforming(fes, mesh)


def test_link_hp_and_mesh_orders():
    mesh = hp2d(3)
    V = H1(mesh, order=4, hpslope=1)
    W = H1(mesh, orderlink=V, relorder=1)
    assert same_space(W, H1(mesh, variableorder=True, relorder=1))
    assert edge_orders(W, mesh) == [p+1 for p in edge_orders(V, mesh)]
    V2 = H1(mesh, variableorder=True, relorder=0)
    W2 = H1(mesh, orderlink=V2, relorder=1)
    assert same_space(W2, W) and consistent(W2) and conforming(W2, mesh)


def test_link_3d_and_low_orders():
    mesh = mesh3d()
    V = H1(mesh, order=2)
    V.SetOrder(mesh.Materials("a"), 4)
    V.Update()
    for rel in [-5, -1, 0, 1]:
        W = H1(mesh, orderlink=V, relorder=rel)
        assert consistent(W) and conforming(W, mesh)
    W = H1(mesh, orderlink=V, relorder=1)
    ref = H1(mesh, order=3); ref.SetOrder(mesh.Materials("a"), 5); ref.Update()
    assert same_space(W, ref)


def test_vectorh1_components_are_linked():
    mesh = mesh2d()
    V = VectorH1(mesh, order=2)
    V.SetOrder(mesh.Materials("b"), 4)
    V.SetOrder(ElementId(VOL, 0), 5)
    V.Update()
    ref = H1(mesh, order=2); ref.SetOrder(mesh.Materials("b"), 4); ref.SetOrder(ElementId(VOL, 0), 5); ref.Update()
    assert V.ndof == 2*ref.ndof and V.GetOrder(ElementId(VOL, 0)) == 5
    assert all(same_space(c, ref) for c in V.components)
    with pytest.raises(Exception, match="orders belong to VectorH1"):
        V.components[1].SetOrder(mesh.Materials("b"), 3)
    with pytest.raises(Exception, match="orders belong to VectorH1"):
        V.components[0].SetOrder(mesh.Materials("b"), 3)
    # a vector space linked to a scalar one, and the other way round
    S = H1(mesh, order=3)
    VV = VectorH1(mesh, orderlink=S, relorder=-1)
    P = H1(mesh, orderlink=V, relorder=-1)
    S.SetOrder(mesh.Materials("a"), 5)
    S.Update(); VV.Update()
    r2 = H1(mesh, order=2); r2.SetOrder(mesh.Materials("a"), 4); r2.Update()
    assert VV.ndof == 2*r2.ndof
    with pytest.raises(Exception, match="orders belong to"):
        VV.SetOrder(mesh.Materials("a"), 3)
    r3 = H1(mesh, order=1); r3.SetOrder(mesh.Materials("b"), 3); r3.SetOrder(ElementId(VOL, 0), 4); r3.Update()
    assert same_space(P, r3)
    # Taylor-Hood with region orders
    X = V * P
    u, p = X.TrialFunction()
    v, q = X.TestFunction()
    a = BilinearForm(InnerProduct(grad(u), grad(v))*dx + div(u)*q*dx + div(v)*p*dx).Assemble()
    assert a.mat.height == V.ndof + P.ndof


def test_link_errors_and_lifetime():
    mesh = mesh2d()
    V = H1(mesh, order=3)
    with pytest.raises(Exception, match="different mesh"):
        H1(mesh2d(), orderlink=V)
    with pytest.raises(Exception, match="orderlink is not available"):
        HDiv(mesh, orderlink=V)
    V.SetOrder(mesh.Materials("b"), 5); V.Update()
    W = H1(mesh, orderlink=V, relorder=-1)
    n = W.ndof
    del V
    gc.collect()
    W.Update()
    assert W.ndof == n and consistent(W)


def test_pickle_linked():
    mesh = mesh2d()
    V = H1(mesh, order=3)
    W = H1(mesh, orderlink=V, relorder=-1)
    V.SetOrder(mesh.Materials("b"), 5)
    V.SetOrder(ElementId(VOL, 0), 4)
    V.Update(); W.Update()
    # together: still linked
    V2, W2 = pickle.loads(pickle.dumps((V, W)))
    assert same_space(V2, V) and same_space(W2, W)
    m2 = V2.mesh
    V2.SetOrder(m2.Materials("a"), 6)
    V2.Update(); W2.Update()
    V.SetOrder(mesh.Materials("a"), 6)
    V.Update(); W.Update()
    assert same_space(V2, V) and same_space(W2, W) and consistent(W2)
    with pytest.raises(Exception, match="orders belong to"):
        W2.SetOrder(m2.Materials("a"), 2)
    # the companion first
    W3, V3 = pickle.loads(pickle.dumps((W, V)))
    V3.SetOrder(V3.mesh.Materials("b"), 2)
    V3.Update(); W3.Update()
    ref = H1(mesh, order=2); ref.SetOrder(mesh.Materials("a"), 5); ref.SetOrder(mesh.Materials("b"), 1)
    ref.SetOrder(ElementId(VOL, 0), 3); ref.Update()
    assert W3.ndof == ref.ndof
    # alone
    W4 = pickle.loads(pickle.dumps(W))
    assert same_space(W4, W)
    W4.Update()
    assert same_space(W4, W)
    V4 = pickle.loads(pickle.dumps(V))
    assert same_space(V4, V)
    # a product space and its GridFunction
    X = V * W
    gf = GridFunction(X)
    gf.components[0].Set(x*y)
    g2 = pickle.loads(pickle.dumps(gf))
    assert len(g2.vec) == len(gf.vec)
    assert abs(Integrate(g2.components[0]**2*dx, g2.space.mesh) - Integrate(gf.components[0]**2*dx, mesh)) < 1e-14
    # vector valued
    VV = VectorH1(mesh, order=2)
    VV.SetOrder(mesh.Materials("b"), 4); VV.Update()
    P = H1(mesh, orderlink=VV, relorder=-1)
    VV2, P2 = pickle.loads(pickle.dumps((VV, P)))
    assert VV2.ndof == VV.ndof and same_space(P2, P)
    VV2.SetOrder(VV2.mesh.Materials("a"), 3)
    VV2.Update(); P2.Update()
    r = H1(mesh, order=1); r.SetOrder(mesh.Materials("b"), 3); r.SetOrder(mesh.Materials("a"), 2); r.Update()
    assert same_space(P2, r)


# -------------------- L2: own orders, linked to H1

def consistent_vol(fes):
    return all(fes.GetFE(el).ndof == len(el.dofs) for el in fes.Elements(VOL))


def ncell(et, p):
    if isinstance(p, int): p = (p, p, p)
    return {ET.TRIG: (p[0]+1)*(p[0]+2)//2, ET.QUAD: (p[0]+1)*(p[1]+1),
            ET.TET: (p[0]+1)*(p[0]+2)*(p[0]+3)//6}[et]


def expected_ndof(mesh, pcell):
    return sum(ncell(el.type, pcell(el)) for el in mesh.Elements(VOL))


@pytest.mark.parametrize("kw", [dict(), dict(all_dofs_together=False), dict(order_policy=ORDER_POLICY.CONSTANT)])
def test_l2_region_and_element_orders(kw):
    for mesh in [mesh2d(), mesh3d()]:
        fes = L2(mesh, order=1, **kw)
        fes.SetOrder(mesh.Materials("b"), 3)
        fes.SetOrder(ElementId(VOL, 0), 4)
        with pytest.raises(Exception, match="order was changed"):
            fes.GetFE(ElementId(VOL, 0))
        fes.Update()
        pcell = lambda el: 4 if el.nr == 0 else 3 if el.mat == "b" else 1
        assert fes.ndof == expected_ndof(mesh, pcell) and consistent_vol(fes)
        assert all(fes.GetOrder(el) == pcell(el) for el in mesh.Elements(VOL))
        mesh.Refine()
        fes.Update()
        assert consistent_vol(fes) and fes.ndof > expected_ndof(mesh, lambda el: 3 if el.mat == "b" else 1)
        # mass matrix works
        u, v = fes.TnT()
        BilinearForm(u*v*dx).Assemble()


def test_l2_direction():
    mesh = MakeStructured2DMesh(quads=True, nx=3, ny=2)
    fes = L2(mesh, order=1)
    fes.SetOrder(mesh.Materials(".*"), (1, 3), direction=(0, 1))
    fes.Update()
    assert fes.ndof == mesh.ne*2*4 and consistent_vol(fes)
    assert sorted(fes.GetOrder(ElementId(VOL, 0))) == [1, 3, 3] or sorted(set(fes.GetOrder(ElementId(VOL, 0)))) == [1, 3]


def test_l2_linked_to_h1():
    mesh = mesh2d()
    V = H1(mesh, order=3)
    Q = L2(mesh, orderlink=V, relorder=-1)
    assert Q.ndof == L2(mesh, order=2).ndof and Q.globalorder == 2
    V.SetOrder(mesh.Materials("b"), 5)
    V.SetOrder(ElementId(VOL, 0), 4)
    V.Update()
    with pytest.raises(Exception, match="order was changed"):
        Q.GetFE(ElementId(VOL, 0))
    Q.Update()
    assert all(Q.GetOrder(el) == V.GetOrder(el)-1 for el in mesh.Elements(VOL))
    assert Q.ndof == expected_ndof(mesh, lambda el: V.GetOrder(el)-1) and consistent_vol(Q)
    with pytest.raises(Exception, match="orders belong to"):
        Q.SetOrder(mesh.Materials("b"), 2)
    with pytest.raises(Exception, match="orders belong to"):
        Q.SetOrder(ElementId(VOL, 1), 2)
    Q.SetOrder(TRIG, 6)          # own layer: an increment of 6 - 2 on triangles
    Q.Update()
    assert Q.ndof == expected_ndof(mesh, lambda el: V.GetOrder(el)-1 + 4) and consistent_vol(Q)


@pytest.mark.parametrize("first", ["creator", "companion"])
def test_l2_link_refinement(first):
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.3))
    V = H1(mesh, order=2)
    Q = L2(mesh, orderlink=V, relorder=-1)
    for el in mesh.Elements(VOL):
        V.SetOrder(el, 2 + el.nr % 3)
    V.Update(); Q.Update()
    for el in mesh.Elements(VOL):
        mesh.SetRefinementFlag(el, el.nr % 3 == 0)
    mesh.Refine()
    for fes in ((V, Q) if first == "creator" else (Q, V)):
        fes.Update()
    assert all(Q.GetOrder(el) == V.GetOrder(el)-1 for el in mesh.Elements(VOL))
    assert Q.ndof == expected_ndof(mesh, lambda el: V.GetOrder(el)-1) and consistent_vol(Q)
    assert len({Q.GetOrder(el) for el in mesh.Elements(VOL)}) == 3


def test_h1_linked_to_l2():
    mesh = mesh2d()
    Q = L2(mesh, order=1)
    V = H1(mesh, orderlink=Q, relorder=1)
    Q.SetOrder(mesh.Materials("b"), 3)
    Q.Update(); V.Update()
    ref = H1(mesh, order=2); ref.SetOrder(mesh.Materials("b"), 4); ref.Update()
    assert V.ndof == ref.ndof
    with pytest.raises(Exception, match="orders belong to"):
        V.SetOrder(mesh.Materials("b"), 3)


def test_l2_hp():
    L = 3
    mesh = hp2d(L)
    ref = L2(mesh, relorder=-1)              # mesh orders - 1
    own = L2(mesh, order=L, hpslope=1)
    V = H1(mesh, order=L+1, hpslope=1)
    Q = L2(mesh, orderlink=V, relorder=-1)
    for fes in (own, Q):
        assert fes.ndof == ref.ndof and consistent_vol(fes)
        assert [tuple(el.dofs) for el in fes.Elements(VOL)] == [tuple(el.dofs) for el in ref.Elements(VOL)]
    assert L2(mesh, order=3, hpslope=0.5).ndof > own.ndof


def test_vector_spaces():
    mesh = mesh2d()
    W = VectorL2(mesh, order=1)
    W.SetOrder(mesh.Materials("b"), 3)
    W.Update()
    s = L2(mesh, order=1); s.SetOrder(mesh.Materials("b"), 3); s.Update()
    assert W.ndof == 2*s.ndof and W.GetOrder(ElementId(VOL, 0)) in (1, 3)
    assert all(W.GetFE(el).ndof == len(el.dofs) for el in W.Elements(VOL))
    with pytest.raises(Exception, match="orders belong to VectorL2"):
        W.components[0].SetOrder(mesh.Materials("b"), 2)

    # Stokes pair with orders on a region
    V = VectorH1(mesh, order=2)
    Q = L2(mesh, orderlink=V, relorder=-2)
    WL = VectorL2(mesh, orderlink=V, relorder=-1)
    V.SetOrder(mesh.Materials("b"), 4)
    for fes in (V, Q, WL): fes.Update()
    assert Q.ndof == expected_ndof(mesh, lambda el: 2 if el.mat == "b" else 0)
    assert WL.ndof == 2*expected_ndof(mesh, lambda el: 3 if el.mat == "b" else 1)
    X = V * Q
    (u, p), (v, q) = X.TnT()
    a = BilinearForm(InnerProduct(grad(u), grad(v))*dx + div(u)*q*dx + div(v)*p*dx).Assemble()
    assert a.mat.height == V.ndof + Q.ndof
    with pytest.raises(Exception, match="orderlink is not available"):
        HDiv(mesh, orderlink=V)


def test_pickle_h1_l2():
    mesh = mesh2d()
    V = H1(mesh, order=3)
    Q = L2(mesh, orderlink=V, relorder=-1)
    V.SetOrder(mesh.Materials("b"), 5)
    V.Update(); Q.Update()
    V2, Q2 = pickle.loads(pickle.dumps((V, Q)))
    assert V2.ndof == V.ndof and Q2.ndof == Q.ndof
    V2.SetOrder(V2.mesh.Materials("a"), 4)
    V2.Update(); Q2.Update()
    assert Q2.ndof == expected_ndof(mesh, lambda el: 4 if el.mat == "b" else 3)
    Q3 = pickle.loads(pickle.dumps(Q))
    assert Q3.ndof == Q.ndof
    gf = GridFunction(Q); gf.Set(x*y)
    g2 = pickle.loads(pickle.dumps(gf))
    assert len(g2.vec) == len(gf.vec)
    fes = L2(mesh, order=1); fes.SetOrder(mesh.Materials("b"), 3); fes.SetOrder(ElementId(VOL, 0), 2); fes.Update()
    f2 = pickle.loads(pickle.dumps(fes))
    assert f2.ndof == fes.ndof
    W = VectorL2(mesh, order=1); W.SetOrder(mesh.Materials("b"), 3); W.Update()
    assert pickle.loads(pickle.dumps(W)).ndof == W.ndof


def test_l2_node_orders_kept():
    mesh = mesh2d()
    fes = L2(mesh, order=1)
    for el in mesh.Elements(VOL):
        if el.nr % 2: fes.SetOrder(el.faces[0], 3)
    fes.SetOrder(mesh.Materials("b"), 2)
    fes.Update()
    pcell = lambda el: 3 if el.nr % 2 else 2 if el.mat == "b" else 1
    assert fes.ndof == expected_ndof(mesh, pcell) and consistent_vol(fes)


# -------------------- HCurl

def tangential_jumps(fes):
    """squared tangential jumps of a function with arbitrary coefficients over the interior facets"""
    u, v = fes.TnT()
    n = specialcf.normal(fes.mesh.dim)
    def tang(w): return w - (w*n)*n
    a = BilinearForm(fes, nonassemble=True)
    a += tang(u-u.Other())*tang(v-v.Other())*dx(skeleton=True)
    gf = GridFunction(fes)
    gf.vec.FV().NumPy()[:] = [((7*i) % 11)/11 - 0.4 for i in range(fes.ndof)]
    res = gf.vec.CreateVector()
    a.Apply(gf.vec, res)
    return abs(InnerProduct(res, gf.vec))


def hcurl_ok(fes):
    return consistent(fes) and tangential_jumps(fes) < 1e-20


def min_rule_edges(mesh, pcell, offset=0):
    pe = {}
    for el in mesh.Elements(VOL):
        for e in el.edges: pe[e.nr] = min(pe.get(e.nr, 99), pcell(el))
    return {i: max(0, p+offset) for i, p in pe.items()}


def test_jump_control_hcurl():
    mesh = mesh2d()
    assert tangential_jumps(HCurl(mesh, order=2)) < 1e-20
    assert tangential_jumps(HCurl(mesh, order=2, discontinuous=True)) > 1e-3


@pytest.mark.parametrize("type1", [False, True])
@pytest.mark.parametrize("make", [mesh2d, mesh3d])
def test_hcurl_region_orders(make, type1):
    mesh = make()
    fes = HCurl(mesh, order=1, type1=type1)
    assert fes.GetOrder(ElementId(VOL, 0)) == 1
    # a region order for all cells is the global order
    fes.SetOrder(mesh.Materials(".*"), 3)
    with pytest.raises(Exception, match="order was changed"):
        fes.GetFE(ElementId(VOL, 0))
    fes.Update()
    assert same_space(fes, HCurl(mesh, order=3, type1=type1))

    fes = HCurl(mesh, order=1, type1=type1)
    fes.SetOrder(mesh.Materials("b"), 3)
    fes.SetOrder(ElementId(VOL, 0), 2)
    fes.Update()
    pcell = lambda el: 2 if el.nr == 0 else 3 if el.mat == "b" else 1
    pe = min_rule_edges(mesh, pcell, -1 if type1 else 0)
    assert all(fes.GetOrder(NodeId(EDGE, i)) == p for i, p in pe.items())
    assert all(fes.GetOrder(el) == pcell(el) for el in mesh.Elements(VOL))
    assert hcurl_ok(fes)
    assert HCurl(mesh, order=1, type1=type1).ndof < fes.ndof < HCurl(mesh, order=3, type1=type1).ndof
    # boundary region, kept by refinement
    fes.SetOrder(mesh.Boundaries("bot"), 4)
    fes.Update()
    bot = {e.nr for el in mesh.Elements(BND) if el.mat == "bot" for e in el.edges}
    assert all(fes.GetOrder(NodeId(EDGE, i)) == (4 - type1) for i in bot)
    assert hcurl_ok(fes)
    n = fes.ndof
    mesh.Refine()
    fes.Update()
    assert fes.ndof > n and hcurl_ok(fes)


def test_hcurl_element_orders_inherited():
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.3))
    fes = HCurl(mesh, order=1)
    pattern = lambda el: 1 + (el.nr * 7) % 3
    for el in mesh.Elements(VOL):
        fes.SetOrder(el, pattern(el))
    fes.Update()
    pe = min_rule_edges(mesh, pattern)
    assert all(fes.GetOrder(NodeId(EDGE, i)) == p for i, p in pe.items()) and hcurl_ok(fes)
    for el in mesh.Elements(VOL):
        mesh.SetRefinementFlag(el, el.nr % 3 == 0)
    mesh.Refine()
    fes.Update()
    assert len({fes.GetOrder(el) for el in mesh.Elements(VOL)}) == 3
    pe = min_rule_edges(mesh, lambda el: fes.GetOrder(el))
    assert all(fes.GetOrder(NodeId(EDGE, i)) == p for i, p in pe.items()) and hcurl_ok(fes)


@pytest.mark.parametrize("make,direction", [(lambda: MakeStructured2DMesh(quads=True, nx=3, ny=3), (0, 1)),
                                            (lambda: MakeStructured3DMesh(hexes=True, nx=2, ny=2, nz=2), (0, 0, 1))])
def test_hcurl_direction(make, direction):
    mesh = make()
    fes = HCurl(mesh, order=1)
    fes.SetOrder(mesh.Materials(".*"), (1, 3), direction=direction)
    fes.Update()
    assert [fes.GetOrder(NodeId(EDGE, e.nr)) for e in mesh.edges] == [3 if along(mesh, e, direction) else 1 for e in mesh.edges]
    assert hcurl_ok(fes)
    assert HCurl(mesh, order=1).ndof < fes.ndof < HCurl(mesh, order=3).ndof


@pytest.mark.parametrize("make", [mesh2d, mesh3d])
def test_hcurl_linked_to_h1(make):
    mesh = make()
    V = H1(mesh, order=3)
    E = HCurl(mesh, orderlink=V, relorder=-1)
    assert same_space(E, HCurl(mesh, order=2)) and E.globalorder == 2
    V.SetOrder(mesh.Materials("b"), 5)
    V.SetOrder(ElementId(VOL, 0), 2)
    V.SetOrder(mesh.Boundaries("bot"), 4)
    V.Update()
    with pytest.raises(Exception, match="order was changed"):
        E.GetFE(ElementId(VOL, 0))
    E.Update()
    used = sorted({e.nr for el in mesh.Elements(VOL) for e in el.edges})
    assert all(E.GetOrder(NodeId(EDGE, i)) == V.GetOrder(NodeId(EDGE, i))-1 for i in used)
    assert all(E.GetOrder(el) == V.GetOrder(el)-1 for el in mesh.Elements(VOL))
    assert hcurl_ok(E)
    with pytest.raises(Exception, match="orders belong to"):
        E.SetOrder(mesh.Materials("b"), 2)
    # gradients of the H1 space are in the HCurl space, also across the changes of the order
    gfv = GridFunction(V)
    gfv.vec.FV().NumPy()[:] = [((5*i) % 13)/13 - 0.3 for i in range(V.ndof)]
    gfe = GridFunction(E)
    gfe.Set(grad(gfv))
    err = Integrate((gfe-grad(gfv))*(gfe-grad(gfv))*dx, mesh)
    assert err < 1e-18 * Integrate(grad(gfv)*grad(gfv)*dx, mesh)
    # pickled together they stay linked
    V2, E2 = pickle.loads(pickle.dumps((V, E)))
    assert same_space(E2, E)
    V2.SetOrder(V2.mesh.Materials("a"), 4)
    V2.Update(); E2.Update()
    V.SetOrder(mesh.Materials("a"), 4)
    V.Update(); E.Update()
    assert same_space(E2, E) and hcurl_ok(E2)


def test_hcurl_hp_and_mesh_orders():
    L = 3
    mesh = hp2d(L)
    ref = HCurl(mesh, variableorder=True, relorder=-1)
    fes = HCurl(mesh, order=L, hpslope=1)
    assert same_space(fes, ref) and hcurl_ok(fes)
    V = H1(mesh, order=L+1, hpslope=1)
    E = HCurl(mesh, orderlink=V, relorder=-1)
    assert same_space(E, ref)


def test_hcurl_node_orders_with_other_calls():
    for mesh in [Mesh(unit_square.GenerateMesh(maxh=0.4)), Mesh(unit_cube.GenerateMesh(maxh=0.7))]:
        et = TRIG if mesh.dim == 2 else TET
        def nodes(fes):
            for e in mesh.edges:
                if e.nr % 2: fes.SetOrder(e, 3)
        a = HCurl(mesh, order=1); a.SetOrder(et, 2); nodes(a); a.Update()
        b = HCurl(mesh, order=1); nodes(b); b.Update(); b.SetOrder(et, 2); b.Update()
        assert same_space(a, b) and hcurl_ok(a)
        assert all(a.GetOrder(NodeId(EDGE, i)) == 3 for i in range(1, mesh.nedge, 2))
        c = HCurl(mesh, order=1); nodes(c)
        with pytest.raises(Exception, match="order was changed"):
            c.GetFE(ElementId(VOL, 0))
        c.UpdateDofTables()
        assert consistent(c)
        f2 = pickle.loads(pickle.dumps(a))
        assert same_space(f2, a)


# -------------------- gradient space of HCurl

def gradient_is_exact(fes):
    """the discrete gradient maps an H1 function to its gradient, and curl grad = 0"""
    mesh = fes.mesh
    gradmat, fesh1 = fes.CreateGradient()
    assert gradmat.height == fes.ndof and gradmat.width == fesh1.ndof
    gfv = GridFunction(fesh1)
    gfv.vec.FV().NumPy()[:] = [((5*i) % 13)/13 - 0.3 for i in range(fesh1.ndof)]
    gfe = GridFunction(fes)
    gfe.vec.data = gradmat * gfv.vec
    err = Integrate((gfe-grad(gfv))*(gfe-grad(gfv))*dx, mesh)
    nrm = Integrate(grad(gfv)*grad(gfv)*dx, mesh)
    u, v = fes.TnT()
    a = BilinearForm(curl(u)*curl(v)*dx).Assemble()
    w = gfe.vec.CreateVector()
    w.data = a.mat * gfe.vec
    return err < 1e-18*nrm and Norm(w) < 1e-9*Norm(gfe.vec)


def setups():
    def regions(f, m): f.SetOrder(m.Materials("b"), 3); f.SetOrder(ElementId(VOL, 0), 2)
    def bnd(f, m): f.SetOrder(m.Boundaries("bot"), 4)
    def nodes(f, m):
        for e in m.edges:
            if e.nr % 2: f.SetOrder(e, 3)
    yield "constant", dict(order=2), None
    yield "constant type1", dict(order=2, type1=True), None
    yield "regions", dict(order=1), regions
    yield "regions type1", dict(order=2, type1=True), regions
    yield "boundary region", dict(order=1), bnd
    yield "node orders", dict(order=1), nodes
    yield "orderinner", dict(order=2, orderinner=3), None
    yield "orderedge", dict(order=3, orderedge=1), None


@pytest.mark.parametrize("make", [mesh2d, mesh3d])
@pytest.mark.parametrize("name", [s[0] for s in setups()])
def test_gradient_space(make, name):
    kw, setup = next((s[1], s[2]) for s in setups() if s[0] == name)
    mesh = make()
    fes = HCurl(mesh, **kw)
    if setup:
        setup(fes, mesh)
        fes.Update()
    assert gradient_is_exact(fes)


@pytest.mark.parametrize("make,L", [(hp2d, 3), (hp3d, 2)])
@pytest.mark.parametrize("slope", [0.5, 1])
def test_gradient_space_hp(make, L, slope):
    mesh = make(L)
    assert gradient_is_exact(HCurl(mesh, order=L+1, hpslope=slope))


def test_gradient_space_direction():
    mesh = MakeStructured3DMesh(hexes=True, nx=2, ny=2, nz=2)
    fes = HCurl(mesh, order=1)
    fes.SetOrder(mesh.Materials(".*"), (1, 3), direction=(0, 0, 1))
    fes.Update()
    assert gradient_is_exact(fes)
