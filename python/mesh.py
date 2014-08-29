import sys
sys.path.append("/opt/netgen/lib")

from libmesh.meshing import *
from libcsg.csg import *


geo = CSGeometry("shaft.geo")
geo.ntlo


param = MeshingParameters()
m1 = GenerateMesh (geo, param)

els = [ i for i in m1.Elements3D() ]
for i in els:
    print (i.vertices)

m1.Save("pymesh.vol")

# mesh = Mesh()
# mesh.Load ("shaft.vol.gz")
# els = mesh.Elements3D()


#cnt = 0
# for el in mesh.Elements3D():
#    print ("el ", cnt, " has vertices " , el.vertices)
#    cnt = cnt+1


