from netgen.meshing import *
from netgen.csg import *

from ngsolve import ngsglobals
ngsglobals.msg_level = 2

geo1 = CSGeometry()
geo1.Add (OrthoBrick( Pnt(0,0,0), Pnt(1,1,1) ))
m1 = geo1.GenerateMesh (maxh=0.1)
m1.Refine()
m1.Refine()


geo2 = CSGeometry()
geo2.Add (Sphere (Pnt(0.5,0.5,0.5), 0.1))
m2 = geo2.GenerateMesh (maxh=0.05)
m2.Refine()
m2.Refine()

print ("********************")
print ("** start merging  **")
print ("********************")

mesh = Mesh()
fd_outside = mesh.Add (FaceDescriptor(surfnr=1,domin=1))
fd_inside = mesh.Add (FaceDescriptor(surfnr=2,domin=2,domout=1))


pmap1 = { }
for e in m1.Elements2D():
    for v in e.vertices:
        if (v not in pmap1):
            pmap1[v] = mesh.Add (m1[v])

for e in m1.Elements2D():
    mesh.Add (Element2D (fd_outside, [pmap1[v] for v in e.vertices]))



pmap2 = { }
for e in m2.Elements2D():
    for v in e.vertices:
        if (v not in pmap2):
            pmap2[v] = mesh.Add (m2[v])

for e in m2.Elements2D():
    mesh.Add (Element2D (fd_inside, [pmap2[v] for v in e.vertices]))


print ("******************")
print ("** merging done **")
print ("******************")


mesh.GenerateVolumeMesh()
mesh.Save ("newmesh.vol")

    
