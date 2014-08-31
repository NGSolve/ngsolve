import sys
sys.path.append("/opt/netgen/lib")

from libmesh.meshing import *
from libcsg.csg import *


geo = CSGeometry("cube.geo")
geo.ntlo


param = MeshingParameters()
# param.maxh = 100
print (param)

m1 = GenerateMesh (geo, param)


for el in m1.Elements3D():
    vi = el.vertices
    for j in vi:
        print (j.nr, m1[j].p)
    print ()


print ("num points = ", len (m1.Points()))

for p in m1.Points():
    print (p.p)


m2 = Mesh()

for p in m1.Points():
    l = p.p
    print (l)
    m2.Add ( MeshPoint (Point(l[0],l[1],l[2])) )
    
print ("Mesh2 is ", m2)



# els = [ i for i in m1.Elements3D() ]
# for i in els:
#     print (i.vertices)

# m1.Save("pymesh.vol")

# mesh = Mesh()
# mesh.Load ("shaft.vol.gz")
# els = mesh.Elements3D()


#cnt = 0
# for el in mesh.Elements3D():
#    print ("el ", cnt, " has vertices " , el.vertices)
#    cnt = cnt+1


