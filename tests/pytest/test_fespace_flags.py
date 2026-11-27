import pickle
import pytest
from ngsolve import *
from ngsolve.comp import GlobalSpace
from netgen.occ import Box, Pnt, OCCGeometry, Glue, X, Y

def make_mesh():
    inner = Box(Pnt(0,0,0), Pnt(0.5,1,1)).mat("inner")
    outer = Box(Pnt(0.5,0,0), Pnt(1,1,1)).mat("outer")
    shape = Glue([inner, outer])
    shape.faces.Min(X).name = "left"
    shape.faces.Max(X).name = "right"
    shape.faces.Min(Y).name = "bot"
    return Mesh(OCCGeometry(shape).GenerateMesh(maxh=0.4))

mesh = make_mesh()

def dofs(fes):
    return fes.ndof, list(fes.FreeDofs())

def test_dirichlet_forms():
    ref = dofs(H1(mesh, order=2, dirichlet="left|bot"))
    assert dofs(H1(mesh, order=2, dirichlet=mesh.Boundaries("left|bot"))) == ref
    assert dofs(H1(mesh, order=2, dirichlet=BND("left|bot"))) == ref
    nrs = [i+1 for i, name in enumerate(mesh.GetBoundaries()) if name in ("left", "bot")]
    assert dofs(H1(mesh, order=2, dirichlet=nrs)) == ref
    fes = H1(mesh, order=2, dirichlet="left|bot")
    assert dofs(pickle.loads(pickle.dumps(fes))) == ref
    # component dirichlet of VectorH1 combines with dirichlet
    v = VectorH1(mesh, order=1, dirichlet="bot", dirichletx="left")
    vref = VectorH1(mesh, order=1, dirichletx="left|bot", dirichlety="bot", dirichletz="bot")
    assert dofs(v) == dofs(vref)

def test_dirichlet_bbnd_forms():
    assert dofs(H1(mesh, order=2, dirichlet_bbnd=".*")) == \
        dofs(H1(mesh, order=2, dirichlet_bbnd=mesh.BBoundaries(".*")))

def test_definedon_forms():
    ref = dofs(H1(mesh, order=2, definedon=mesh.Materials("inner")))
    assert dofs(H1(mesh, order=2, definedon="inner")) == ref
    assert dofs(H1(mesh, order=2, definedon=VOL("inner"))) == ref
    assert dofs(H1(mesh, order=2, definedon=["inner"])) == ref
    nr = list(mesh.GetMaterials()).index("inner") + 1
    assert dofs(H1(mesh, order=2, definedon=[nr])) == ref
    fes = H1(mesh, order=2, definedon="inner")
    assert dofs(pickle.loads(pickle.dumps(fes))) == ref

def test_definedon_dict_shim():
    ref = dofs(H1(mesh, order=2, definedon=mesh.Materials("inner") | mesh.Boundaries("right")))
    assert dofs(H1(mesh, order=2, definedon={VOL: mesh.Materials("inner"), BND: mesh.Boundaries("right")})) == ref
    assert dofs(H1(mesh, order=2, definedon={VOL: "inner", BND: "right"})) == ref
    fes = H1(mesh, order=2, definedon={VOL: mesh.Materials("inner"), BND: mesh.Boundaries("right")})
    assert dofs(pickle.loads(pickle.dumps(fes))) == ref
    # only boundary: no volume dofs
    assert dofs(H1(mesh, order=2, definedon={BND: mesh.Boundaries("right")})) == \
        dofs(H1(mesh, order=2, definedon=mesh.Boundaries("right")))
    with pytest.raises(Exception, match="keys must be"):
        H1(mesh, definedon={"inner": 1})
    with pytest.raises(Exception, match="needs a VOL region"):
        H1(mesh, definedon={VOL: mesh.Boundaries("right")})

def test_order_policy_flag():
    # H1 derives its policy from the orders, L2 keeps the flag
    for op in [ORDER_POLICY.CONSTANT, ORDER_POLICY.NODETYPE, ORDER_POLICY.VARIABLE]:
        assert L2(mesh, order=2, order_policy=op).order_policy == op

def test_orderlink_flag():
    V = H1(mesh, order=3)
    W = H1(mesh, orderlink=V, relorder=-1)
    assert W.ndof == H1(mesh, order=2).ndof
    other = Mesh(unit_cube.GenerateMesh(maxh=0.5))
    with pytest.raises(Exception, match="different mesh"):
        H1(other, orderlink=V)

def test_globalspace_basis():
    fes = GlobalSpace(mesh, order=1, basis=CF((1, x, y)))
    assert fes.ndof == 3
    gf = GridFunction(H1(mesh, order=1))
    gf.Set(x)
    assert GlobalSpace(mesh, order=1, basis=gf).ndof == 1

def test_bfi_definedon():
    fes = H1(mesh, order=1)
    a1 = BilinearForm(fes)
    a1 += BFI("mass", coef=1, definedon=[list(mesh.GetMaterials()).index("inner")+1])
    a1.Assemble()
    a2 = BilinearForm(fes)
    a2 += fes.TrialFunction()*fes.TestFunction()*dx("inner")
    a2.Assemble()
    gf = GridFunction(fes); gf.Set(1)
    assert abs(InnerProduct(a1.mat*gf.vec, gf.vec) - 0.5) < 1e-10
    assert abs(InnerProduct(a2.mat*gf.vec, gf.vec) - 0.5) < 1e-10
