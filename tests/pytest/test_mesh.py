import operator
import pytest
from ngsolve import *
from netgen.csg import *
ngsglobals.msg_level = 0
from netgen.geom2d import *

def test_multiple_meshes_refine():
    mesh = Mesh(unit_cube.GenerateMesh(maxh=1))
    mesh.Refine()
    parents = [mesh.GetParentVertices(v) for v in range(mesh.nv)]
    print(parents)

    mesh2 = Mesh(unit_cube.GenerateMesh(maxh=0.5))
    mesh2.Refine()

    for v in range(mesh.nv):
        assert parents[v] == mesh.GetParentVertices(v)

def test_call_operator_with_vectorization():
    mesh = Mesh(unit_cube.GenerateMesh(maxh=1))
    p = mesh(0.5,0.5,0.5)
    p2 = mesh([0.5, 0.1],0.5,0.5)

def test_neighbours():
    geo = CSGeometry()

    bot = Plane(Pnt(0,0,0), Vec(0,0,-1)).bc("bot")
    top = Plane(Pnt(10,10,10), Vec(0,0,1)).bc("top")
    box = OrthoBrick(Pnt(0,0,-1), Pnt(10,10,11)).bc("rest") * top * bot

    plate = OrthoBrick(Pnt(3,3,4.8), Pnt(7,7,5.2))
    cutx = Plane(Pnt(5,5,5), Vec(1,0,0)).bc("cutx")
    cuty = Plane(Pnt(5,5,5), Vec(0,1,0)).bc("cuty")

    cond1 = plate * cutx * cuty
    cond2 = plate * cutx - cuty
    cond3 = plate * cuty - cutx
    cond4 = plate - cutx - cuty

    geo.Add(cond1.mat("cond1"))
    geo.Add(cond2.mat("cond2"))
    geo.Add(cond3.mat("cond3"))
    geo.Add(cond4.mat("cond4"))

    geo.Add((box-plate).mat("air"))

    mesh = Mesh(geo.GenerateMesh(maxh=2))

    xpl = mesh.Boundaries("cutx")
    nmats = [mesh.GetMaterials()[i] for i, val in enumerate(xpl.Neighbours(VOL).Mask()) if val]
    print(nmats)
    assert nmats == [ "cond1", "cond2", "cond3", "cond4"]
    cond12 = mesh.Materials("cond1")
    print(set(mesh.GetBoundaries()[i] for i, val in enumerate(cond12.Boundaries().Mask()) if val))
    assert set((mesh.GetBoundaries()[i] for i, val in enumerate(cond12.Boundaries().Mask()) if val)) == {"default", "cuty", "cutx"}
    print(set(mesh.GetBoundaries()[i] for i, val in enumerate(mesh.Materials("air").Boundaries().Mask()) if val))
    assert not any([bnd in set(mesh.GetBoundaries()[i] for i, val in enumerate(mesh.Materials("air").Boundaries().Mask()) if val) for bnd in ["cuty", "cutx"]])

    import numpy
    midedge = numpy.prod([reg.Neighbours(BBND) for reg in mesh.Materials("cond.*").Split()])
    print(midedge.Mask())
    # CSG has one edge region per adjacent face pair, all on the mid edge x=y=5
    assert sum(midedge.Mask()) >= 1
    for el in mesh.Elements(BBND):
        if midedge.Mask()[el.index]:
            assert all(abs(mesh[v].point[0]-5) < 1e-8 and abs(mesh[v].point[1]-5) < 1e-8 for v in el.vertices)

    midedge2 = mesh.Boundaries("cuty").Boundaries() * mesh.Boundaries("cutx").Boundaries()
    print(midedge2.Mask())
    assert midedge == midedge2
    print(mesh.Materials("cond1").Neighbours(VOL).Mask())
    assert mesh.Materials("cond1").Neighbours(VOL) == mesh.Materials("cond2|cond3|air")

def test_neighbours2d():
    geo = CSG2d()
    top = Rectangle((0.2,0.6), (0.8, 0.8), bc="outer", bottom="default", mat="top")
    base = Rectangle((0,0), (1, 0.6), bc="outer")
    chip = Solid2d([(0.5,0.15), (0.65,0.3), (0.5,0.45), (0.35,0.3)], mat="chip")
    geo.Add((base-chip).Mat("base"))
    geo.Add(top)
    geo.Add(chip)
    mesh = Mesh(geo.GenerateMesh())
    Draw(mesh)
    print(mesh.Materials("base").Boundaries().Mask())
    print(mesh.Materials("chip").Boundaries().Mask())
    print(mesh.Materials("top").Boundaries().Mask())
    assert mesh.Materials("base").Boundaries() * mesh.Materials("top").Boundaries() == mesh.Boundaries("default")
    assert mesh.Materials("base").Boundaries() * mesh.Materials("chip").Boundaries() == mesh.Boundaries("")

def test_region_operators():
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.3))
    r = mesh.Boundaries("left") + mesh.Boundaries("right")
    assert r.VB() == BND and sum(r.Mask()) == 2
    assert sum((r - mesh.Boundaries("left")).Mask()) == 1
    assert sum((r * mesh.Boundaries("left")).Mask()) == 1
    assert sum((r + "top").Mask()) == 3

@pytest.mark.parametrize("op", [operator.add, operator.sub, operator.mul])
def test_region_operators_mixed_vb(op):
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.3))
    with pytest.raises(Exception, match="cannot combine VOL and BND"):
        op(mesh.Materials(".*"), mesh.Boundaries("left"))

def test_region_operators_other_mesh():
    m1 = Mesh(unit_square.GenerateMesh(maxh=0.3))
    m2 = Mesh(unit_square.GenerateMesh(maxh=0.3))
    with pytest.raises(Exception, match="different meshes"):
        m1.Boundaries("left") + m2.Boundaries("left")

def test_neighbours_all_codims_occ():
    import netgen.occ as occ
    box = occ.Box((0,0,0), (1,1,1))
    box.faces.name = "side"
    box.faces.Min(occ.Z).name = "bot"
    box.faces.Max(occ.Z).name = "top"
    box.edges.Min(occ.Z).name = "botedge"
    box.vertices.Min(occ.X+occ.Y+occ.Z).name = "corner"
    box.solids.name = "air"
    ball = occ.Sphere((0.5,0.5,0.5), 0.25)
    ball.faces.name = "sphere"
    ball.solids.name = "ball"
    mesh = Mesh(occ.OCCGeometry(occ.Glue([box - ball, ball])).GenerateMesh(maxh=0.3))
    assert mesh.Materials("ball").Boundaries() == mesh.Boundaries("sphere")
    assert mesh.Boundaries("sphere").Neighbours(VOL) == mesh.Materials("air|ball")
    edge_faces = mesh.BBoundaries("botedge").Neighbours(BND)
    assert sum(edge_faces.Mask()) == 2 and edge_faces * mesh.Boundaries("bot") == mesh.Boundaries("bot")
    assert mesh.BBoundaries("botedge").Neighbours(VOL) == mesh.Materials("air")
    assert mesh.BBBoundaries("corner").Neighbours(VOL) == mesh.Materials("air")
    corner_faces = mesh.BBBoundaries("corner").Neighbours(BND)
    assert sum(corner_faces.Mask()) == 3 and corner_faces * mesh.Boundaries("bot") == mesh.Boundaries("bot")
    assert sum(mesh.BBBoundaries("corner").Neighbours(BBND).Mask()) == 3
    assert mesh.Materials("ball").Neighbours(BBBND) * mesh.BBBoundaries("corner") == mesh.BBBoundaries("")

def test_neighbours_edge_contact_occ():
    import netgen.occ as occ
    b1 = occ.Box((0,0,0), (1,1,1)); b1.solids.name = "b1"
    b2 = occ.Box((1,1,0), (2,2,1)); b2.solids.name = "b2"
    mesh = Mesh(occ.OCCGeometry(occ.Glue([b1, b2])).GenerateMesh(maxh=0.5))
    shared = mesh.Materials("b1").Neighbours(BBND) * mesh.Materials("b2").Neighbours(BBND)
    assert sum(shared.Mask()) == 1
    assert sum(mesh.Materials("b1").Neighbours(VOL).Mask()) == 0

def test_neighbours_2d_occ():
    import netgen.occ as occ
    rect = occ.Rectangle(2, 1).Face()
    rect.edges.name = "outer"
    rect.vertices.Min(occ.X+occ.Y).name = "pin"
    rect.faces.name = "outer"
    circ = occ.Circle((1, 0.5), 0.2).Face()
    circ.edges.name = "circ"
    circ.faces.name = "inner"
    mesh = Mesh(occ.OCCGeometry(occ.Glue([rect - circ, circ]), dim=2).GenerateMesh(maxh=0.2))
    assert mesh.Materials("inner").Boundaries() == mesh.Boundaries("circ")
    assert mesh.Materials("outer").Boundaries() == mesh.Boundaries("outer|circ")
    assert mesh.BBoundaries("pin").Neighbours(VOL) == mesh.Materials("outer")
    assert sum(mesh.BBoundaries("pin").Neighbours(BND).Mask()) == 2

def _box_ball_mesh():
    import netgen.occ as occ
    box = occ.Box((0,0,0), (1,1,1))
    box.faces.name = "side"
    box.faces.Min(occ.Z).name = "bot"
    box.vertices.Min(occ.X+occ.Y+occ.Z).name = "corner"
    box.solids.name = "air"
    ball = occ.Sphere((0.5,0.5,0.5), 0.25)
    ball.faces.name = "sphere"
    ball.solids.name = "ball"
    return Mesh(occ.OCCGeometry(occ.Glue([box - ball, ball])).GenerateMesh(maxh=0.3))

def test_definedon_closed():
    # a space is defined on the closure of its definedon regions
    mesh = _box_ball_mesh()
    for reg in [mesh.Materials("ball"), mesh.Boundaries("bot"), mesh.BBoundaries(".*"),
                mesh.Materials("ball") | mesh.BBBoundaries("corner"),
                {VOL: mesh.Materials("ball"), BND: mesh.Boundaries("bot")}]:
        fes = H1(mesh, order=2, definedon=reg)
        if isinstance(reg, dict):
            reg = reg[VOL] | reg[BND]
        if isinstance(reg, Region):
            reg = MixedRegion(reg)
        closed = reg.Closure()
        for vb in [VOL, BND, BBND, BBBND]:
            for el in mesh.Elements(vb):
                assert (len(fes.GetDofNrs(ElementId(vb, el.nr))) > 0) == closed[vb].Mask()[el.index]

def test_fespace_definedon_property():
    mesh = _box_ball_mesh()
    assert H1(mesh, order=2).definedon == ~MixedRegion(mesh)
    fes = H1(mesh, order=2, definedon=mesh.Materials("ball"))
    assert fes.definedon == MixedRegion(mesh.Materials("ball")).Closure()
    for vb in [VOL, BND, BBND, BBBND]:
        assert fes.GetDefinedOnRegion(vb) == fes.definedon[vb]

def test_fespace_element_without_region():
    import netgen.meshing as nm
    m = nm.Mesh(3)
    pts = [m.Add(nm.MeshPoint(nm.Pnt(*p))) for p in [(0,0,0),(1,0,0),(0,1,0),(0,0,1)]]
    m.Add(nm.Element3D(1, pts))
    fd = m.Add(nm.FaceDescriptor(surfnr=1, domin=1, bc=1))
    for t in [(0,1,2),(0,1,3),(0,2,3),(1,2,3)]:
        m.Add(nm.Element2D(fd, [pts[i] for i in t]))
    m.Add(nm.Element1D([pts[0], pts[1]], index=0))
    m.SetMaterial(1, "a")
    mesh = Mesh(m)
    assert [e.index for e in mesh.Elements(BBND)] == [-1]
    fes = H1(mesh, order=2)
    assert fes.ndof == 10
    assert len(fes.GetDofNrs(ElementId(BBND, 0))) == 0
    fes = H1(mesh, order=2, definedon=mesh.Materials("a"))
    assert len(fes.GetDofNrs(ElementId(BBND, 0))) == 0

def test_dirichlet_element_without_region():
    import netgen.meshing as nm
    m = nm.Mesh(3)
    pts = [m.Add(nm.MeshPoint(nm.Pnt(*p))) for p in [(0,0,0),(1,0,0),(0,1,0),(0,0,1)]]
    m.Add(nm.Element3D(1, pts))
    fd = m.Add(nm.FaceDescriptor(surfnr=1, domin=1, bc=1))
    for t in [(0,1,2),(0,1,3),(0,2,3),(1,2,3)]:
        m.Add(nm.Element2D(fd, [pts[i] for i in t]))
    m.Add(nm.Element1D([pts[0], pts[1]], index=0))
    m.Add(nm.Element1D([pts[2], pts[3]], index=1))
    mesh = Mesh(m)
    assert [e.index for e in mesh.Elements(BBND)] == [-1, 0]
    fes = H1(mesh, order=1, dirichlet_bbnd=".*")
    assert list(fes.FreeDofs()) == [True, True, False, False]
    fes = H1(mesh, order=1, dirichlet=mesh.BBoundaries(".*"))
    assert list(fes.FreeDofs()) == [True, True, False, False]

def test_mixedregion_setops():
    import pickle
    mesh = _box_ball_mesh()
    clamp = mesh.Boundaries("bot") | mesh.BBBoundaries("corner")
    assert isinstance(clamp, MixedRegion)
    assert clamp[BND] == mesh.Boundaries("bot")
    assert clamp[BBBND] == mesh.BBBoundaries("corner")
    assert clamp[VOL] == mesh.Materials("")
    assert clamp == mesh.MixedRegion(BND="bot", BBBND="corner")
    assert hash(clamp) == hash(mesh.BBBoundaries("corner") | mesh.Boundaries("bot"))
    assert (clamp & mesh.Boundaries("bot|side")) == MixedRegion(mesh.Boundaries("bot"))
    assert (clamp - mesh.Boundaries("bot")) == MixedRegion(mesh.BBBoundaries("corner"))
    assert mesh.MixedRegion().Empty() and not clamp.Empty()
    assert ~mesh.MixedRegion() == mesh.MixedRegion(VOL=".*", BND=".*", BBND=".*", BBBND=".*")
    assert isinstance(mesh.Boundaries("bot") + mesh.Boundaries("side"), Region)
    m2, c2 = pickle.loads(pickle.dumps((mesh, clamp)))
    assert c2 == m2.Boundaries("bot") | m2.BBBoundaries("corner")
    other = _box_ball_mesh()
    with pytest.raises(Exception, match="different meshes"):
        clamp | other.Boundaries("bot")

def test_mixedregion_topology():
    mesh = _box_ball_mesh()
    ball = MixedRegion(mesh.Materials("ball"))
    assert ball.Boundary() == MixedRegion(mesh.Boundaries("sphere"))
    cl = ball.Closure()
    assert cl[VOL] == mesh.Materials("ball") and cl[BND] == mesh.Boundaries("sphere")
    assert cl.IsClosed() and not ball.IsClosed()
    air = MixedRegion(mesh.Materials("air")).Closure()
    assert air[BND] == mesh.Materials("air").Boundaries()
    assert air[BBBND] == mesh.Materials("air").Neighbours(BBBND)
    star = MixedRegion(mesh.BBBoundaries("corner")).Star()
    assert star[VOL] == mesh.Materials("air")
    assert sum(star[BND].Mask()) == 3 and sum(star[BBND].Mask()) == 3
    assert star[BBBND] == mesh.BBBoundaries("corner")
    # interface of the two volumes
    iface = MixedRegion(mesh.Materials("air")).Closure() & MixedRegion(mesh.Materials("ball")).Closure()
    assert iface[VOL] == mesh.Materials("") and iface[BND] == mesh.Boundaries("sphere")

def test_mixedregion_fespace_dirichlet_property():
    mesh = _box_ball_mesh()
    reg = mesh.Boundaries("bot") | mesh.BBBoundaries("corner")
    fes = H1(mesh, order=2, dirichlet=reg)
    assert fes.dirichlet == reg
    d = fes.dirichlet
    d = d - MixedRegion(mesh.Boundaries("bot"))
    assert d != reg and fes.dirichlet == reg
    with pytest.raises(AttributeError):
        fes.dirichlet = d

@pytest.mark.parametrize("space", [H1, HCurl])
def test_mixedregion_fespace_flags(space):
    mesh = _box_ball_mesh()
    def dofs(**kwargs):
        fes = space(mesh, order=2, **kwargs)
        return fes.ndof, list(fes.FreeDofs())
    assert dofs(dirichlet=MixedRegion(mesh.Boundaries("bot"))) == dofs(dirichlet="bot")
    assert dofs(dirichlet=mesh.Boundaries("bot") | mesh.BBBoundaries("corner")) == \
        dofs(dirichlet="bot", dirichlet_bbbnd="corner")
    assert dofs(dirichlet=MixedRegion(mesh.BBBoundaries("corner"))) == dofs(dirichlet_bbbnd="corner")
    assert dofs(dirichlet=mesh.BBBoundaries("corner")) == dofs(dirichlet_bbbnd="corner")
    assert dofs(dirichlet=mesh.BBoundaries(".*")) == dofs(dirichlet_bbnd=".*")
    assert dofs(dirichlet=mesh.BBoundaries(".*"), dirichlet_bbbnd="corner") == \
        dofs(dirichlet_bbnd=".*", dirichlet_bbbnd="corner")
    with pytest.raises(Exception, match="VOL regions are not supported"):
        space(mesh, order=2, dirichlet=mesh.Materials("ball"))
    assert dofs(dirichlet=MixedRegion(mesh.Boundaries("bot")), dirichlet_bbbnd="corner") == \
        dofs(dirichlet="bot", dirichlet_bbbnd="corner")
    assert dofs(definedon=MixedRegion(mesh.Materials("ball"))) == dofs(definedon=mesh.Materials("ball"))
    assert dofs(definedon=MixedRegion(mesh.Boundaries("sphere"))) == dofs(definedon=mesh.Boundaries("sphere"))
    assert dofs(definedon=mesh.Materials("ball") | mesh.Boundaries("bot")) == \
        dofs(definedon={VOL: mesh.Materials("ball"), BND: mesh.Boundaries("bot")})
    with pytest.raises(Exception, match="VOL regions are not supported"):
        space(mesh, dirichlet=MixedRegion(mesh.Materials("ball")))

def test_mixedregion_getdofs():
    mesh = _box_ball_mesh()
    mr = mesh.Materials("ball") | mesh.Boundaries("bot") | mesh.BBBoundaries("corner")
    for space in (H1, HCurl, L2):
        fes = space(mesh, order=2)
        union = fes.GetDofs(mr[VOL]) | fes.GetDofs(mr[BND]) | fes.GetDofs(mr[BBBND])
        assert list(fes.GetDofs(mr)) == list(union)
        assert list(fes.GetDofs(MixedRegion(mesh.Boundaries("bot")))) == list(fes.GetDofs(mesh.Boundaries("bot")))
    u1, u2 = (H1(mesh, order=1) * H1(mesh, order=2)).TrialFunction()
    assert list(u2.GetDofs(mr)) == list(u2.GetDofs(mr[VOL]) | u2.GetDofs(mr[BND]) | u2.GetDofs(mr[BBBND]))

def test_mixedregion_indicator():
    import pickle
    mesh = _box_ball_mesh()
    ind = (mesh.Materials("ball") | mesh.Boundaries("bot")).IndicatorCF()
    assert Integrate(ind, mesh, VOL) == pytest.approx(Integrate(1, mesh, definedon=mesh.Materials("ball")))
    assert Integrate(ind, mesh, BND) == pytest.approx(1)
    assert Integrate(ind.Compile(), mesh, VOL) == pytest.approx(Integrate(ind, mesh, VOL))
    assert ind(mesh(0.5, 0.5, 0.5)) == 1 and ind(mesh(0.1, 0.1, 0.9)) == 0
    assert Integrate((ind*x).Diff(x), mesh, VOL) == pytest.approx(Integrate(ind, mesh, VOL))
    mesh2, ind2 = pickle.loads(pickle.dumps((mesh, ind)))
    assert ind2(mesh2(0.5, 0.5, 0.5)) == 1 and ind2(mesh2(0.1, 0.1, 0.9)) == 0

def test_mixed_regioncf():
    import pickle
    mesh = _box_ball_mesh()
    vball = Integrate(1, mesh, definedon=mesh.Materials("ball"))
    cf = mesh.RegionCF({mesh.Materials("ball"): 2, mesh.Materials("air"): 1,
                        mesh.Boundaries("bot"): 5}, default=0)
    assert Integrate(cf, mesh, VOL) == pytest.approx(1 + vball)
    assert Integrate(cf, mesh, BND) == pytest.approx(5)
    assert Integrate(cf.Compile(), mesh, BND) == pytest.approx(5)
    assert pickle.loads(pickle.dumps(cf))(mesh(0.5, 0.5, 0.5)) == 2
    cf2 = mesh.RegionCF({mesh.Materials("ball"): x*x, mesh.Boundaries("bot"): y}, default=0)
    assert Integrate(cf2.Diff(x), mesh, VOL) == pytest.approx(Integrate(2*x, mesh, definedon=mesh.Materials("ball")))
    assert Integrate(cf2.Diff(y), mesh, BND) == pytest.approx(1)
    cm = mesh.RegionCF({mesh.Materials("ball") | mesh.Boundaries("sphere"): 3}, default=0)
    assert Integrate(cm, mesh, BND) == pytest.approx(3*Integrate(1, mesh, definedon=mesh.Boundaries("sphere")))

    # proxies in the values: same matrix as the per-VorB form
    fes = H1(mesh, order=2)
    u, v = fes.TnT()
    cfu = mesh.RegionCF({mesh.Materials(".*"): u, mesh.Boundaries("bot"): 3*u})
    a1 = BilinearForm(cfu*v*dx + cfu*v*ds(definedon=mesh.Boundaries("bot"))).Assemble()
    a2 = BilinearForm(u*v*dx + 3*u*v*ds("bot")).Assemble()
    gf = GridFunction(fes)
    gf.Set(x*y+z)
    r = gf.vec.CreateVector()
    r.data = a1.mat*gf.vec - a2.mat*gf.vec
    assert Norm(r) < 1e-12

    with pytest.raises(Exception, match="needs Region or MixedRegion"):
        mesh.RegionCF({"ball": 1})
    with pytest.raises(Exception, match="evaluated on BND"):
        Integrate(mesh.RegionCF({mesh.Materials("ball"): 1}), mesh, BND)
    with pytest.raises(Exception, match="same dimension"):
        mesh.RegionCF({mesh.Materials("ball"): 1, mesh.Boundaries("bot"): CF((1,2))})
    # same geometry, same region counts: only the mesh check catches it
    other = _box_ball_mesh()
    with pytest.raises(Exception, match="another mesh"):
        mesh.RegionCF({other.Materials("ball"): 1})
    with pytest.raises(Exception, match="another mesh"):
        mesh.RegionCF({other.Materials("ball") | other.Boundaries("bot"): 1})
    with pytest.raises(Exception, match="another mesh"):
        mesh.RegionCF(VOL, {other.Materials("ball"): 1})

def _brute_neighbours(mesh, vb1, vb2):
    # regions of vb2 touching each region of vb1, from the vertex sets of the elements
    from collections import defaultdict
    hi, lo = (vb1, vb2) if int(vb1) < int(vb2) else (vb2, vb1)
    his = [(frozenset(v.nr for v in el.vertices), el.index) for el in mesh.Elements(hi)]
    v2hi = defaultdict(set)
    for i, (vs, _) in enumerate(his):
        for v in vs: v2hi[v].add(i)
    res = defaultdict(set)
    for el in mesh.Elements(lo):
        vs = frozenset(v.nr for v in el.vertices)
        for i in set.intersection(*(v2hi[v] for v in vs)):
            if vb1 == hi: res[his[i][1]].add(el.index)
            else: res[el.index].add(his[i][1])
    return res

def _check_neighbours(mesh):
    vbs = [VOL, BND, BBND, BBBND][:mesh.dim+1]
    for vb1 in vbs:
        for vb2 in vbs:
            if vb1 == vb2: continue
            exp = _brute_neighbours(mesh, vb1, vb2)
            for i, reg in enumerate(mesh.Region(vb1, ".*").Split()):
                got = set(j for j, b in enumerate(reg.Neighbours(vb2).Mask()) if b)
                assert got == exp.get(i, set()), (vb1, vb2, i)

def test_neighbours_curved_pickled_loaded(tmp_path):
    import pickle
    mesh = _box_ball_mesh()
    mesh.Curve(3)
    _check_neighbours(mesh)
    assert H1(mesh, dirichlet=MixedRegion(mesh.Boundaries("bot")).Closure()).ndof > 0
    _check_neighbours(pickle.loads(pickle.dumps(_box_ball_mesh())))
    _box_ball_mesh().ngmesh.Save(str(tmp_path / "boxball.vol"))
    _check_neighbours(Mesh(str(tmp_path / "boxball.vol")))
    ngmesh = _box_ball_mesh().ngmesh
    ngmesh.SecondOrder()
    _check_neighbours(Mesh(ngmesh))

def test_neighbours_follow_mesh_changes():
    ngmesh = _box_ball_mesh().ngmesh
    mesh = Mesh(ngmesh)
    _check_neighbours(mesh)
    # move the elements of the first boundary region to the second one
    for el in ngmesh.Elements2D():
        if el.index == 1: el.index = 2
    fresh = Mesh(ngmesh)
    assert not fresh.Materials(".*").Neighbours(BND).Mask()[0]
    _check_neighbours(fresh)
    mesh.Refine()
    assert not mesh.Materials(".*").Neighbours(BND).Mask()[0]
    _check_neighbours(mesh)

def test_neighbours_quads_hexes_prisms():
    from ngsolve.meshes import MakeStructured2DMesh, MakeStructured3DMesh
    _check_neighbours(MakeStructured2DMesh(quads=True, nx=3, ny=3))
    _check_neighbours(MakeStructured3DMesh(hexes=True, nx=2, ny=2, nz=2))
    _check_neighbours(MakeStructured3DMesh(hexes=False, prism=True, nx=2, ny=2, nz=2))

def test_region_python_operators():
    mesh = _box_ball_mesh()
    bot, side = mesh.Boundaries("bot"), mesh.Boundaries("side")
    mr = bot | mesh.BBBoundaries("corner")
    # same VorB: stays a Region
    assert isinstance(bot | side, Region) and (bot | side) == bot + side
    assert isinstance(mr, MixedRegion) and isinstance(bot | mr, MixedRegion) and isinstance(mr | bot, MixedRegion)
    assert not (mr == 3) and mr != None and mr in [1, mr] and not (bot == 3)
    # Region and MixedRegion never compare equal, like their hashes
    assert MixedRegion(bot) != bot and bot != MixedRegion(bot)
    assert {mr: 1}[mesh.BBBoundaries("corner") | bot] == 1
    with pytest.raises(TypeError):
        bot | 3

def test_mixedregion_misuse():
    mesh, other = _box_ball_mesh(), _box_ball_mesh()
    ind = (mesh.Materials("ball") | mesh.Boundaries("bot")).IndicatorCF()
    with pytest.raises(Exception, match="another mesh"):
        Integrate(ind, other, VOL)
    with pytest.raises(Exception, match="no region given"):
        mesh.RegionCF({MixedRegion(mesh.Materials("typo")): 1}, default=0)

def _box_rail_mesh():
    import netgen.occ as occ
    box = occ.Box((0,0,0), (1,1,1))
    box.faces.name = "side"
    box.faces.Min(occ.Z).name = "bot"
    box.faces.Max(occ.Z).name = "top"
    box.edges.Max(occ.X+occ.Z).name = "rail"
    box.vertices.Max(occ.X+occ.Y+occ.Z).name = "corner"
    box.solids.name = "air"
    ball = occ.Sphere((0.5,0.5,0.5), 0.25)
    ball.faces.name = "sphere"
    ball.solids.name = "ball"
    return Mesh(occ.OCCGeometry(occ.Glue([box - ball, ball])).GenerateMesh(maxh=0.3))

def test_mixedregion_descriptor():
    import pickle
    mesh = _box_rail_mesh()
    d = BND("bot") | BBND("rail") | BBBND("corner")
    assert str(d) == "BND(bot) | BBND(rail) | BBBND(corner)"
    assert mesh[d] == mesh.Boundaries("bot") | mesh.BBoundaries("rail") | mesh.BBBoundaries("corner")
    assert mesh[d | BND("top")] == mesh[d] | mesh.Boundaries("top")
    assert str(d[BND]) == "BND(bot)" and d[VOL] is None
    assert str(pickle.loads(pickle.dumps(d))) == str(d)
    # same VorB: stays a RegionDescriptor
    same = BND("bot") | BND("top")
    assert not hasattr(same, "__getitem__") and mesh[same] == mesh.Boundaries("bot|top")

    def dofs(fes): return fes.ndof, list(fes.FreeDofs())
    ref = H1(mesh, order=2, dirichlet="bot", dirichlet_bbnd="rail", dirichlet_bbbnd="corner")
    fes = H1(mesh, order=2, dirichlet=d)
    assert dofs(fes) == dofs(ref)
    assert dofs(pickle.loads(pickle.dumps(fes))) == dofs(ref)
    assert dofs(H1(mesh, order=2, dirichlet=BND("bot"))) == dofs(H1(mesh, order=2, dirichlet="bot"))
    assert dofs(H1(mesh, order=2, definedon=VOL("ball") | BND("bot"))) == \
        dofs(H1(mesh, order=2, definedon={VOL: mesh.Materials("ball"), BND: mesh.Boundaries("bot")}))
    assert dofs(H1(mesh, order=2, definedon=VOL("ball"))) == dofs(H1(mesh, order=2, definedon="ball"))
    with pytest.raises(Exception, match="VOL regions are not supported"):
        H1(mesh, dirichlet=VOL("ball") | BND("bot"))

    gf = GridFunction(fes)
    gf[d] = 3
    assert gf(mesh(0.3, 0.3, 0)) == pytest.approx(3) and gf(mesh(1, 0.5, 1)) == pytest.approx(3)
    assert gf(mesh(1, 1, 1)) == pytest.approx(3) and gf(mesh(0.1, 0.1, 0.9)) == pytest.approx(0)

    cf = mesh.RegionCF({VOL("ball"): 2, VOL("air"): 1, BND("bot") | BBND("rail"): 5}, default=0)
    assert Integrate(cf, mesh, BND) == pytest.approx(5) and Integrate(cf, mesh, BBND) == pytest.approx(5)
    assert Integrate(cf, mesh, VOL) == pytest.approx(1 + Integrate(1, mesh, definedon=mesh.Materials("ball")))

def test_mixedregion_descriptor_solve():
    mesh = _box_rail_mesh()
    d = BND("bot") | BBND("rail") | BBBND("corner")
    fref = H1(mesh, order=2, dirichlet="bot", dirichlet_bbnd="rail", dirichlet_bbbnd="corner")
    ur, vr = fref.TnT()
    a = BilinearForm(grad(ur)*grad(vr)*dx + ur*vr*ds("top")).Assemble()
    f = LinearForm(vr*dx).Assemble()
    gref = GridFunction(fref)
    gref.vec.data = a.mat.Inverse(fref.FreeDofs()) * f.vec

    fes = H1(mesh, order=2)
    u, v = fes.TnT()
    lhs, rhs = grad(u)*grad(v)*dx + u*v*ds("top"), 1*v*dx
    def err(gf): return sqrt(Integrate((gf-gref)**2, mesh))
    assert err(Solve(lhs == rhs, u[d] == 0)) < 1e-12
    for pre in (preconditioners.Local, preconditioners.BDDC, preconditioners.MultiGrid):
        assert err(Solve(lhs == rhs, u[d] == 0, pre(), solvers.CGSolver.Creator(maxiter=2000, tol=1e-12))) < 1e-10
    gin = Solve(lhs == rhs, u[d] == 2)
    assert gin(mesh(0.3, 0.3, 0)) == pytest.approx(2) and gin(mesh(1, 1, 1)) == pytest.approx(2)

    bc = u[d] == 0
    assert str(bc.vbn) == str(d)
    with pytest.raises(Exception, match="one is needed"):
        bc.vb
    single = u[BND("bot")] == 0
    assert single.vb == BND and single.name == "bot" and str(single.vbn) == "BND(bot)"

if __name__ == "__main__":
    test_neighbours2d()
    test_neighbours()
