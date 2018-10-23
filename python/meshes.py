from netgen.meshing import *
from netgen.csg import *
import ngsolve

def MakeQuadMesh(nx=10, ny=10, mapping = lambda x,y : (x,y) ):
    mesh = Mesh()
    mesh.dim=2

    pids = []
    for i in range(ny+1):
        for j in range(nx+1):
            x,y = mapping (j/nx, i/ny)
            pids.append (mesh.Add (MeshPoint(Pnt(x,y,0))))

    mesh.Add (FaceDescriptor(surfnr=1,domin=1,bc=1))
            
    for i in range(ny):
        for j in range(nx):
            base = i * (nx+1) + j
            pnum = [base,base+1,base+nx+2,base+nx+1]
            elpids = [pids[p] for p in pnum]
            mesh.Add (Element2D(1,elpids))

    for i in range(nx):
        mesh.Add(Element1D([pids[i], pids[i+1]], index=1))
    for i in range(ny):
        mesh.Add(Element1D([pids[i*(nx+1)+nx], pids[(i+1)*(nx+1)+nx]], index=2))
    for i in range(nx):
        mesh.Add(Element1D([pids[ny*(nx+1)+i+1], pids[ny*(nx+1)+i]], index=3))
    for i in range(ny):
        mesh.Add(Element1D([pids[(i+1)*(nx+1)], pids[i*(nx+1)]], index=4))
    mesh.SetBCName(0, "bottom")        
    mesh.SetBCName(1, "right")        
    mesh.SetBCName(2, "top")        
    mesh.SetBCName(3, "left")        

    return ngsolve.Mesh(mesh)
