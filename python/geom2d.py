from nglib.meshing import *
from nglib.geom2d import *

geom = SplineGeometry()
geom.Load("square.in2d")

param = MeshingParameters()
param.maxh = 0.05

m1 = geom.GenerateMesh (param)