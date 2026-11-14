import pytest
from ngsolve import *
from ngsolve.meshes import MakeStructured2DMesh, MakeStructured3DMesh
from ngsolve.krylovspace import CGSolver

ngsglobals.msg_level = 0


def Mesh2D(**kwargs):
    return MakeStructured2DMesh(quads=False, nx=4, ny=3, **kwargs)


def Mesh3D(**kwargs):
    return MakeStructured3DMesh(hexes=False, nx=3, ny=2, nz=2, **kwargs)


def Rows(table):
    return [list(row) for row in table]


def FreeDofs(fes, coupling=False):
    return set(i for i, b in enumerate(fes.FreeDofs(coupling)) if b)


@pytest.mark.parametrize("mesh", [Mesh2D(), Mesh3D()])
@pytest.mark.parametrize("space", [lambda m: H1(m, order=4, dirichlet="left|bottom"),
                                   lambda m: HCurl(m, order=2, dirichlet="left|bottom"),
                                   lambda m: HDiv(m, order=2, dirichlet="left|bottom"),
                                   lambda m: VectorH1(m, order=2, dirichlet="left") * H1(m, order=1),
                                   lambda m: L2(m, order=2) * FacetFESpace(m, order=2, dirichlet="left")])
def test_named_blocks(mesh, space):
    fes = space(mesh)
    free = FreeDofs(fes)
    dim = mesh.dim

    # the open nodes are a partition of the free dofs
    rows = Rows(fes.CreateSmoothingBlocks(["vertex", "edge", "face", "cell"]))
    alldofs = [d for row in rows for d in row]
    assert len(alldofs) == len(set(alldofs))
    assert set(alldofs) == free

    nnodes = {"vertex": mesh.nv, "edge": mesh.nedge, "face": mesh.nface}
    nnodes["facet"] = nnodes["edge"] if dim == 2 else nnodes["face"]
    nnodes["cell"] = mesh.ne if dim == 3 else 0
    nnodes["inner"] = mesh.ne
    nnodes["element"] = mesh.ne
    for name, anchor in (("vertex", "vertex"), ("edge", "edge"), ("face", "face"), ("cell", "cell"),
                         ("facet", "facet"), ("inner", "inner"), ("element", "element"), ("vertexedge", "vertex"),
                         ("vertexpatch", "vertex"), ("edgepatch", "edge"), ("facepatch", "face")):
        rows = Rows(fes.CreateSmoothingBlocks(name))
        assert len(rows) == nnodes[anchor]
        for row in rows:
            assert row == sorted(set(row))
            assert set(row) <= free
        assert Rows(fes.CreateSmoothingBlocks(blocktype=[name])) == rows

    # a patch contains the smaller blocks around the same node
    v = Rows(fes.CreateSmoothingBlocks("vertex"))
    ve = Rows(fes.CreateSmoothingBlocks("vertexedge"))
    vp = Rows(fes.CreateSmoothingBlocks("vertexpatch"))
    for a, b, c in zip(v, ve, vp):
        assert set(a) <= set(b) <= set(c)
    assert set(d for row in vp for d in row) == free

    # codimension names
    facet = Rows(fes.CreateSmoothingBlocks("facet"))
    inner = Rows(fes.CreateSmoothingBlocks("inner"))
    assert facet == Rows(fes.CreateSmoothingBlocks("edge" if dim == 2 else "face"))
    assert sorted(inner) == sorted(Rows(fes.CreateSmoothingBlocks("face" if dim == 2 else "cell")))
    assert Rows(fes.CreateSmoothingBlocks(["inner", "face", "cell"])) == \
        inner + Rows(fes.CreateSmoothingBlocks("cell" if dim == 2 else "face"))
    for el in fes.Elements():
        assert set(inner[el.nr]) <= set(el.dofs)

    # element blocks hold the element dofs
    rows = Rows(fes.CreateSmoothingBlocks("element"))
    for el in fes.Elements():
        assert set(rows[el.nr]) == set(d for d in el.dofs if d >= 0) & free

    # several types: rows of the single types, in the order of the list
    e = Rows(fes.CreateSmoothingBlocks("edge"))
    assert Rows(fes.CreateSmoothingBlocks(["vertexpatch", "edge"])) == vp + e
    assert Rows(fes.CreateSmoothingBlocks(["edge", "vertexpatch", "edge"])) == e + vp


def test_errors():
    fes = H1(Mesh2D(), order=2)
    with pytest.raises(Exception):
        fes.CreateSmoothingBlocks("vertpatch")
    with pytest.raises(Exception):
        fes.CreateSmoothingBlocks(["vertex", "bogus"])
    with pytest.raises(Exception):
        fes.CreateSmoothingBlocks("element:periodic")
    with pytest.raises(Exception):
        fes.CreateSmoothingBlocks("vertexpatch:a,b")


@pytest.mark.parametrize("mesh", [Mesh2D(), Mesh3D()])
def test_condense(mesh):
    fes = H1(mesh, order=5, dirichlet="left")
    for kwargs in ({"blocktype": "vertexpatch"}, {}, {"blocktype": 3}):
        full = set(d for row in fes.CreateSmoothingBlocks(**kwargs) for d in row)
        cond = set(d for row in fes.CreateSmoothingBlocks(condense=True, **kwargs) for d in row)
        assert full == FreeDofs(fes)
        assert cond == FreeDofs(fes, True)
        assert Rows(fes.CreateSmoothingBlocks(eliminate_internal=True, **kwargs)) == \
            Rows(fes.CreateSmoothingBlocks(condense=True, **kwargs))


def test_h1_default_blocks():
    # default of H1: vertex with its edges, faces, inner dofs
    fes = H1(Mesh3D(), order=4, dirichlet="left")
    default = Rows(fes.CreateSmoothingBlocks())
    named = Rows(fes.CreateSmoothingBlocks(["vertexedge", "face", "inner"]))
    assert Rows(fes.CreateSmoothingBlocks(blocktype=4)) == default
    assert default == named
    nv, nf = fes.mesh.nv, fes.mesh.nface
    assert default[:nv] == Rows(fes.CreateSmoothingBlocks("vertexedge"))
    assert default[nv:nv+nf] == Rows(fes.CreateSmoothingBlocks("face"))
    # 2D: no face blocks, inner blocks numbered like the elements
    fes = H1(Mesh2D(), order=4, dirichlet="left")
    assert Rows(fes.CreateSmoothingBlocks()) == Rows(fes.CreateSmoothingBlocks(["vertexedge", "inner"]))


def test_threads():
    mesh = Mesh3D()
    for fes in (H1(mesh, order=3), HCurl(mesh, order=2), HDiv(mesh, order=1) * L2(mesh, order=1)):
        for blocktype in ("vertexpatch", "edgepatch", ["vertexedge", "face", "inner"], "element", None):
            sequential = Rows(fes.CreateSmoothingBlocks(blocktype))
            with TaskManager():
                parallel = Rows(fes.CreateSmoothingBlocks(blocktype))
            assert [sorted(row) for row in sequential] == [sorted(row) for row in parallel]
            if blocktype is not None:
                assert sequential == parallel


@pytest.mark.parametrize("mesh", [Mesh2D(periodic_x=True, periodic_y=True),
                                  MakeStructured3DMesh(hexes=False, nx=3, periodic_x=True, periodic_y=True)])
def test_periodic(mesh):
    fes = Periodic(H1(mesh, order=3))
    free = FreeDofs(fes)
    for name in ("vertex", "edge", "vertexpatch", "edgepatch"):
        plain = Rows(fes.CreateSmoothingBlocks(name))
        merged = Rows(fes.CreateSmoothingBlocks(name + ":periodic"))
        assert len(merged) < len(plain)
        assert set(d for row in merged for d in row) == set(d for row in plain for d in row)
        assert all(row == sorted(set(row)) for row in merged)
    # one vertex dof per merged vertex block
    vertexdofs = set(d for row in fes.CreateSmoothingBlocks("vertex") for d in row)
    merged = Rows(fes.CreateSmoothingBlocks("vertex:periodic"))
    assert len(merged) == len(vertexdofs)
    assert sorted(merged) == sorted([d] for d in vertexdofs)
    assert set(d for row in fes.CreateSmoothingBlocks("vertexpatch:periodic") for d in row) == free


def test_dofset():
    mesh = Mesh2D()
    fes = HDiv(mesh, order=2)
    lo = set(d for row in fes.CreateSmoothingBlocks("vertexpatch:hdivlo") for d in row)
    assert lo == set(range(mesh.nedge))
    rows = Rows(fes.CreateSmoothingBlocks("edge:hdivlo"))
    assert rows == [[i] for i in range(mesh.nedge)]


@pytest.mark.parametrize("gs", [False, True])
def test_local_preconditioner(gs):
    mesh = Mesh3D()
    fes = H1(mesh, order=3, dirichlet="left|bottom")
    u, v = fes.TnT()

    def Iterations(**kwargs):
        a = BilinearForm(grad(u)*grad(v)*dx)
        f = LinearForm(v*dx)
        pre = Preconditioner(a, "local", GS=gs, **kwargs)
        a.Assemble()
        f.Assemble()
        inv = CGSolver(a.mat, pre.mat, maxiter=500, tol=1e-10, printrates=False)
        gfu = GridFunction(fes)
        gfu.vec.data = inv * f.vec
        return inv.iterations

    point = Iterations()
    patch = Iterations(blocktype="vertexpatch")
    assert patch < point
    assert Iterations(block=True) == Iterations(blocktype=4)
    assert Iterations(blocktype=["vertexedge", "face", "inner"]) == Iterations(blocktype=4)
    assert Iterations(blockcreator=lambda fes: fes.CreateSmoothingBlocks("vertexpatch")) == patch
