import numpy as np
import pytest
from ngsolve import *
from netgen.meshing import Mesh as NGMesh, MeshPoint, Element3D, Element2D, FaceDescriptor, Pnt

# proxy.forms: how each factor of a proxy value transforms from the reference element

V, E, F, L, C = FORM.VERTEX, FORM.EDGE, FORM.FACET, FORM.ELEMENT, FORM.COMP


def single_element_mesh(dim, A=None):
    A = np.eye(dim) if A is None else A
    m = NGMesh(dim=dim)
    if dim == 3:
        pis = [m.Add(MeshPoint(Pnt(*(A @ np.array(p))))) for p in [(0,0,0),(1,0,0),(0,1,0),(0,0,1)]]
        m.SetMaterial(1, "dom")
        m.Add(Element3D(1, pis))
        m.Add(FaceDescriptor(surfnr=1, domin=1, bc=1))
        for f in [(1,2,3),(0,3,2),(0,1,3),(0,2,1)]:
            m.Add(Element2D(1, [pis[i] for i in f]))
    else:
        pis = [m.Add(MeshPoint(Pnt(*(A @ np.array(p)), 0))) for p in [(0,0),(1,0),(0,1)]]
        m.SetMaterial(1, "dom")
        m.Add(Element2D(1, pis))
    return Mesh(m)


def factor_maps(F, forward):
    J = np.linalg.det(F) if F.shape[0] == F.shape[1] else np.sqrt(np.linalg.det(F.T @ F))
    Fp = np.linalg.pinv(F)
    if forward:
        return { V: (None, 1), E: (Fp.T, 1), F_: (F/J, 1), L: (None, 1/J), C: (np.eye(F.shape[0]), 1) }
    return { V: (None, 1), E: (F.T, 1), F_: (J*Fp, 1), L: (None, J), C: (np.eye(F.shape[0]), 1) }

F_ = F


def transform(forms, detpower, Fmat, val, forward):
    maps = factor_maps(Fmat, forward)
    J = np.linalg.det(Fmat) if Fmat.shape[0] == Fmat.shape[1] else np.sqrt(np.linalg.det(Fmat.T @ Fmat))
    mats, w = [], J**detpower if forward else J**(-detpower)
    for f in forms:
        M, s = maps[f]
        w *= s
        if M is not None:
            mats.append(M)
    v = np.array(val, dtype=float)
    if mats:
        v = v.reshape([M.shape[1] for M in mats])
        for i, M in enumerate(mats):
            v = np.moveaxis(np.tensordot(M, v, axes=([1], [i])), 0, i)
    return w * np.array(v).flatten()


def jacobian(mesh, mp, vb):
    d = mesh.dim
    de = d if vb == VOL else d-1
    return np.array(specialcf.JacobianMatrix(d, de)(mp)).reshape(d, de)


def mapping_error(mref, mphys, mkfes, op, vb, pts_ref, pts_phys):
    """max deviation of the physical values from the transformation given by proxy.forms"""
    fr, fp = mkfes(mref), mkfes(mphys)
    p = op(fp.TrialFunction())
    forms, detpower = p.forms, p.forms_detpower
    gr, gp = GridFunction(fr), GridFunction(fp)
    err = scale = 0
    for i in range(fr.ndof):
        gr.vec[:] = 0
        gp.vec[:] = 0
        gr.vec[i] = 1
        gp.vec[i] = 1
        for qr, qp in zip(pts_ref, pts_phys):
            mpr, mpp = mref(*qr, vb), mphys(*qp, vb)
            vr = np.atleast_1d(np.array(op(gr)(mpr), dtype=float))
            vp = np.atleast_1d(np.array(op(gp)(mpp), dtype=float))
            vhat = transform(forms, detpower, jacobian(mref, mpr, vb), vr, False)
            ex = transform(forms, detpower, jacobian(mphys, mpp, vb), vhat, True)
            err = max(err, np.abs(ex-vp).max())
            scale = max(scale, np.abs(vp).max())
    return err / scale


Id = lambda u: u
Op = lambda name: (lambda u: u.Operator(name))

cases3d = [
    ("H1",              lambda m: H1(m, order=3), Id, (V,), True),
    ("H1 grad",         lambda m: H1(m, order=3), grad, (E,), True),
    ("H1 hesse",        lambda m: H1(m, order=3), Op("hesse"), (E,E), False),
    ("H1 dim=3",        lambda m: H1(m, order=2, dim=3), Id, (C,V), True),
    ("H1 dim=3 Grad",   lambda m: H1(m, order=2, dim=3), Grad, (C,E), True),
    ("HCurl",           lambda m: HCurl(m, order=2), Id, (E,), True),
    ("HCurl curl",      lambda m: HCurl(m, order=2), curl, (F,), True),
    ("HCurl grad",      lambda m: HCurl(m, order=2), Op("grad"), (E,E), False),
    ("HDiv",            lambda m: HDiv(m, order=2), Id, (F,), True),
    ("HDiv div",        lambda m: HDiv(m, order=2), div, (L,), True),
    ("HDiv grad",       lambda m: HDiv(m, order=2), Op("grad"), (E,F), False),
    ("L2",              lambda m: L2(m, order=2), Id, (V,), True),
    ("L2 dual",         lambda m: L2(m, order=2), Op("dual"), (L,), True),
    ("VectorL2",        lambda m: VectorL2(m, order=2), Id, (C,V), True),
    ("VectorL2 piola",  lambda m: VectorL2(m, order=2, piola=True), Id, (F,), True),
    ("VectorL2 piola div", lambda m: VectorL2(m, order=2, piola=True), div, (L,), True),
    ("VectorL2 covariant", lambda m: VectorL2(m, order=2, covariant=True), Id, (E,), True),
    ("VectorL2 covariant curl", lambda m: VectorL2(m, order=2, covariant=True), curl, (F,), True),
    ("HCurlCurl",       lambda m: HCurlCurl(m, order=1), Id, (E,E), True),
    ("HCurlCurl inc",   lambda m: HCurlCurl(m, order=2), Op("inc"), (F,F), False),
    ("HCurlCurl curl",  lambda m: HCurlCurl(m, order=1), curl, (E,F), False),
    ("HDivDiv",         lambda m: HDivDiv(m, order=1), Id, (F,F), True),
    ("HDivDiv div",     lambda m: HDivDiv(m, order=1), div, (F,L), False),
    ("HCurlDiv",        lambda m: HCurlDiv(m, order=1), Id, (E,F), True),
    ("HCurlDiv div",    lambda m: HCurlDiv(m, order=1), div, (E,L), False),
    ("VectorH1",        lambda m: VectorH1(m, order=2), Id, (C,V), True),
    ("VectorH1 Grad",   lambda m: VectorH1(m, order=2), Grad, (C,E), True),
    ("VectorH1 hesse",  lambda m: VectorH1(m, order=2), Op("hesse"), (C,E,E), False),
    ("VectorValued(HCurl)", lambda m: VectorValued(HCurl(m, order=1)), Id, (C,E), True),
    ("VectorValued(HCurl) curl", lambda m: VectorValued(HCurl(m, order=1)), curl, (C,F), True),
    ("VectorValued(HCurlCurl)", lambda m: VectorValued(HCurlCurl(m, order=1)), Id, (C,E,E), True),
    ("MatrixValued(H1)", lambda m: MatrixValued(H1(m, order=1)), Id, (C,C,V), True),
    ("MatrixValued(H1) Grad", lambda m: MatrixValued(H1(m, order=1)), Grad, (C,C,E), True),
    ("MatrixValued(H1) sym", lambda m: MatrixValued(H1(m, order=1), symmetric=True), Id, (C,C,V), True),
]

cases3d_bnd = [
    ("H1 trace",        lambda m: H1(m, order=2), lambda u: u.Trace(), (V,), True),
    ("HCurl trace",     lambda m: HCurl(m, order=2), lambda u: u.Trace(), (E,), True),
    ("HCurlCurl trace", lambda m: HCurlCurl(m, order=1), lambda u: u.Trace(), (E,E), True),
    ("SurfaceL2",       lambda m: SurfaceL2(m, order=1), Id, (V,), True),
    ("SurfaceL2 dual_mapping", lambda m: SurfaceL2(m, order=1, dual_mapping=True), Id, (L,), True),
    ("TangentialFacet trace", lambda m: TangentialFacetFESpace(m, order=1), lambda u: u.Trace(), (E,), True),
    ("FacetFESpace trace", lambda m: FacetFESpace(m, order=1), lambda u: u.Trace(), (V,), True),
]

cases2d = [
    ("H1 grad",         lambda m: H1(m, order=3), grad, (E,), True),
    ("HCurl",           lambda m: HCurl(m, order=2), Id, (E,), True),
    ("HCurl curl",      lambda m: HCurl(m, order=2), curl, (L,), True),
    ("HDiv",            lambda m: HDiv(m, order=2), Id, (F,), True),
    ("HDiv div",        lambda m: HDiv(m, order=2), div, (L,), True),
    ("HCurlCurl",       lambda m: HCurlCurl(m, order=2), Id, (E,E), True),
    ("HCurlCurl inc",   lambda m: HCurlCurl(m, order=2), Op("inc"), (L,L), False),
    ("HCurlCurl curl",  lambda m: HCurlCurl(m, order=2), curl, (E,L), False),
    ("HDivDiv",         lambda m: HDivDiv(m, order=2), Id, (F,F), True),
    ("HDivDiv div",     lambda m: HDivDiv(m, order=2), div, (F,L), False),
    ("HCurlDiv",        lambda m: HCurlDiv(m, order=2), Id, (E,F), True),
    ("HCurlDiv div",    lambda m: HCurlDiv(m, order=2), div, (E,L), False),
]

A3 = np.array([[1.3, 0.4, -0.2], [0.1, 0.8, 0.5], [-0.3, 0.2, 1.7]])
A2 = np.array([[1.3, 0.4], [-0.2, 0.9]])
pts3 = [np.array([0.2, 0.3, 0.1]), np.array([0.1, 0.15, 0.5])]
pts2 = [np.array([0.2, 0.3]), np.array([0.1, 0.5])]
pts3_bnd = [np.array([0.0, 0.3, 0.2]), np.array([0.0, 0.1, 0.6])]


def ids(cases):
    return [c[0] for c in cases]


@pytest.mark.parametrize("case", cases3d, ids=ids(cases3d))
def test_forms_3d(case):
    name, mkfes, op, forms, exact = case
    mref, mphys = single_element_mesh(3), single_element_mesh(3, A3)
    p = op(mkfes(mphys).TrialFunction())
    assert p.forms == forms
    assert p.forms_exact == exact
    assert mapping_error(mref, mphys, mkfes, op, VOL, pts3, [A3 @ q for q in pts3]) < 1e-12


@pytest.mark.parametrize("case", cases3d_bnd, ids=ids(cases3d_bnd))
def test_forms_3d_boundary(case):
    name, mkfes, op, forms, exact = case
    mref, mphys = single_element_mesh(3), single_element_mesh(3, A3)
    p = op(mkfes(mphys).TrialFunction())
    assert p.forms == forms
    assert p.forms_exact == exact
    assert mapping_error(mref, mphys, mkfes, op, BND, pts3_bnd, [A3 @ q for q in pts3_bnd]) < 1e-12


@pytest.mark.parametrize("case", cases2d, ids=ids(cases2d))
def test_forms_2d(case):
    name, mkfes, op, forms, exact = case
    mref, mphys = single_element_mesh(2), single_element_mesh(2, A2)
    p = op(mkfes(mphys).TrialFunction())
    assert p.forms == forms
    assert p.forms_exact == exact
    assert mapping_error(mref, mphys, mkfes, op, VOL, pts2, [A2 @ q for q in pts2]) < 1e-12


@pytest.mark.parametrize("case", cases3d, ids=ids(cases3d))
def test_forms_exact_on_curved_element(case):
    name, mkfes, op, forms, exact = case
    mref, mcurved = single_element_mesh(3), single_element_mesh(3)
    deform = GridFunction(VectorH1(mcurved, order=2))
    deform.Set((0.1*y*z+0.05*x*x, 0.08*x*z-0.04*y*y, 0.06*x*y+0.03*z*z))
    mcurved.SetDeformation(deform)
    err = mapping_error(mref, mcurved, mkfes, op, VOL, pts3, pts3)
    if exact:
        assert err < 1e-12
    else:
        assert err > 1e-6


dual = Op("dual")
cases_dual = [
    ("H1 dual",         lambda m: H1(m, order=4), (L,), 0),
    ("HCurl dual",      lambda m: HCurl(m, order=3), (F,), 0),
    ("HDiv dual",       lambda m: HDiv(m, order=2), (E,), 0),
    ("L2 dual",         lambda m: L2(m, order=2), (L,), 0),
    ("VectorH1 dual",   lambda m: VectorH1(m, order=4), (C,L), 0),
    ("HCurlCurl dual",  lambda m: HCurlCurl(m, order=2), (F,F), 1),
    ("HDivDiv dual",    lambda m: HDivDiv(m, order=1), (E,E), 1),
    ("HCurlDiv dual",   lambda m: HCurlDiv(m, order=1), (F,E), 1),
]
# duals evaluated on boundary faces, types relative to the face
cases_dual_bnd = [c for c in cases_dual if c[0] in ("H1 dual", "HCurl dual", "HCurlCurl dual")]


@pytest.mark.parametrize("case", cases_dual, ids=ids(cases_dual))
def test_forms_dual(case):
    name, mkfes, forms, detpower = case
    mref, mphys = single_element_mesh(3), single_element_mesh(3, A3)
    p = dual(mkfes(mphys).TrialFunction())
    assert p.forms == forms
    assert p.forms_detpower == detpower
    assert mapping_error(mref, mphys, mkfes, dual, VOL, pts3, [A3 @ q for q in pts3]) < 1e-12


@pytest.mark.parametrize("case", cases_dual_bnd, ids=ids(cases_dual_bnd))
def test_forms_dual_boundary(case):
    name, mkfes, forms, detpower = case
    mref, mphys = single_element_mesh(3), single_element_mesh(3, A3)
    assert mapping_error(mref, mphys, mkfes, dual, BND, pts3_bnd, [A3 @ q for q in pts3_bnd]) < 1e-12


def test_forms_compound_and_untyped():
    mesh = single_element_mesh(3)
    u, p = (HCurl(mesh, order=1)*H1(mesh, order=1)).TrialFunction()
    assert u.forms == (E,) and grad(p).forms == (E,)
    assert NumberSpace(mesh).TrialFunction().forms is None
    assert HDiv(mesh, order=1).TrialFunction().Trace().forms is None    # normal vector times density


def test_form_enum():
    assert [f.degree(3) for f in (V, E, F, L)] == [0, 1, 2, 3]
    assert [f.degree(2) for f in (V, E, F, L)] == [0, 1, 1, 2]
    assert [f.star for f in (V, E, F, L)] == [L, F, E, V]
    assert [f.from_top for f in (V, E, F, L, C)] == [False, False, True, True, False]
    assert [f.has_index for f in (V, E, F, L, C)] == [False, True, True, False, True]
    assert C.degree(3) is None
