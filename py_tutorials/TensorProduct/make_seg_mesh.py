from netgen.meshing import *
from netgen.csg import *

def AddEdgeEls (p1, dx, facenr, n, pids,mesh):
    for i in range(n):
        base = p1 + i*dx
        pnum = [base, base+dx]
        elpids = [pids[p] for p in pnum]
        mesh.Add (Element1D(elpids,index=1))

def SegMesh(n,x0,x1,per):
    mesh = Mesh(dim=1)
    pids = []
    for i in range(n+1):
        pids.append (mesh.Add (MeshPoint(Pnt(x0+(x1-x0)*i/n, 0, 0))))
    AddEdgeEls(0,1,1,n,pids,mesh)
    mesh.Add (Element0D( pids[0], index=1))
    mesh.Add (Element0D( pids[n], index=2))
    if per == 1:
        mesh.AddPointIdentification(pids[0],pids[n],1,2)
    return mesh