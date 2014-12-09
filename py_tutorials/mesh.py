from netgen.meshing import *
from netgen.csg import *

geo = CSGeometry("shaft.geo")

param = MeshingParameters()
param.maxh = 10
print (param)

m1 = GenerateMesh (geo, param)


import exportNeutral
exportNeutral.Export (m1, "shaft.mesh")



Save (m1, "mesh.vol", geo)
