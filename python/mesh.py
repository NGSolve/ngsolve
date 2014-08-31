import sys
sys.path.append("/opt/netgen/lib")

from libmesh.meshing import *
from libcsg.csg import *


geo = CSGeometry("shaft.geo")

param = MeshingParameters()
param.maxh = 10
print (param)

m1 = GenerateMesh (geo, param)


import exportNeutral
exportNeutral.Export (m1, "shaft.mesh")



