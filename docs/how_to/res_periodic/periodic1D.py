from netgen.meshing import *

mesh = Mesh(dim=1)
pids = []
n = 50
for i in range(n+1):
    pids.append (mesh.Add (MeshPoint(Pnt(i/n, 0, 0))))
for i in range(n):
    mesh.Add(Element1D([pids[i],pids[i+1]],index=1))
mesh.Add (Element0D( pids[0], index=1))
mesh.Add (Element0D( pids[n], index=2))
mesh.AddPointIdentification(pids[0],pids[n],1,2)
