# tests for high order H(div) pyramids: normal continuity to tets and hexes on a hybrid mesh,
# polynomial containment, divergence into the L2 space
import pytest
from netgen.meshing import Mesh as NgMesh, MeshPoint, Element3D, Element2D, FaceDescriptor, Pnt
from ngsolve import *

def hybrid_mesh(perm):
    # pyramid with base (0,0,0)-(1,1,0) and apex (0.5,0.5,1); tets on the trig faces, hex below the base
    pts = {0:(0,0,0), 1:(1,0,0), 2:(1,1,0), 3:(0,1,0), 4:(0.5,0.5,1),
           5:(0.5,-1,0.6), 6:(2,0.5,0.6), 7:(0.5,2,0.6), 8:(-1,0.5,0.6),
           9:(0,0,-1), 10:(1,0,-1), 11:(1,1,-1), 12:(0,1,-1)}
    m = NgMesh(dim=3)
    pid = {}
    for k in perm:      # order in which points are added defines the vertex numbering
        pid[k] = m.Add(MeshPoint(Pnt(*pts[k])))
    idx = m.AddRegion("vol", dim=3)
    m.Add(Element3D(idx, [pid[i] for i in [0,1,2,3,4]]))
    for f,q in [((0,1,4),5),((1,2,4),6),((2,3,4),7),((3,0,4),8)]:
        m.Add(Element3D(idx, [pid[f[0]],pid[f[1]],pid[f[2]],pid[q]]))
    m.Add(Element3D(idx, [pid[i] for i in [9,10,11,12,0,1,2,3]]))
    fd = m.Add(FaceDescriptor(surfnr=1, domin=idx, bc=1))
    bnd = [(0,1,5),(1,4,5),(4,0,5), (1,2,6),(2,4,6),(4,1,6), (2,3,7),(3,4,7),(4,2,7), (3,0,8),(0,4,8),(4,3,8),
           (9,12,11,10),(9,10,1,0),(10,11,2,1),(11,12,3,2),(12,9,0,3)]
    for f in bnd:
        m.Add(Element2D(fd, [pid[i] for i in f]))
    return Mesh(m)

@pytest.mark.parametrize("order", [0,1,2,3])
@pytest.mark.parametrize("RT", [False, True])
@pytest.mark.parametrize("perm", [list(range(13)), [4,3,2,1,0,8,7,6,5,12,11,10,9],
                                  [2,0,4,1,3,5,6,7,8,9,10,11,12], [6,11,0,4,9,3,1,12,2,10,7,8,5]])
def test_hdiv_pyramid(order, RT, perm):
    mesh = hybrid_mesh(perm)
    fes = HDiv(mesh, order=order, RT=RT, dgjumps=True)
    u,v = fes.TnT()
    n = specialcf.normal(3)
    intorder = 2*order+6
    # normal continuity across all interior facets
    a = BilinearForm(fes)
    a += (u*n-u.Other()*n)*(v*n-v.Other()*n)*dx(skeleton=True)
    a.Assemble()
    assert a.mat.AsVector().Norm() < 1e-10
    # polynomials of degree 'order' are in the space
    gf = GridFunction(fes)
    cfs = [CF((x**order,0,0)), CF((0,y**order,0)), CF((0,0,z**order))]
    if order >= 1:
        cfs.append(CF((x*y**(order-1), z**order, x**order)))
    for cf in cfs:
        gf.Set(cf)
        assert Integrate((gf-cf)**2, mesh, order=intorder) < 1e-20
    # divergence of the pyramid basis functions lands in the L2 space (order-1, or order for RT)
    # (not checked on the hex: HDiv hexes are always RT-type)
    l2 = L2(mesh, order=order if RT else max(order-1,0))
    gfd = GridFunction(l2)
    nohex = [el.type != ET.HEX for el in mesh.Elements(VOL)]
    pyr = [el for el in fes.Elements(VOL) if el.type == ET.PYRAMID][0]
    for i in pyr.dofs:
        gf.vec[:] = 0; gf.vec[i] = 1
        gfd.Set(div(gf))
        err = Integrate((gfd-div(gf))**2, mesh, order=intorder, element_wise=True)
        assert max(e for e,keep in zip(err,nohex) if keep) < 1e-20
