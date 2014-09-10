from nglib.meshing import *
from nglib.geom2d import *

geom = SplineGeometry()
SplineGeometry.Plot = plotgeom
SplineGeometry.ShowPoints = plotpointindex
SplineGeometry.ShowDomains = plotdomainindex

# Define Points
pi1 = geom.AppendPoint(0,0)
pi2 = geom.AppendPoint(1,0)
pi3 = geom.AppendPoint(1,0.5)
pi4 = geom.AppendPoint(1,1)
pi5 = geom.AppendPoint(0.5,1)
pi6 = geom.AppendPoint(0,1)

# Define Segments
geom.AppendSegment(pi1,pi2)
geom.AppendSegment(pi2,pi3)
geom.AppendSegment(pi3,pi4,pi5)
geom.AppendSegment(pi5,pi6)
geom.AppendSegment(pi6,pi1)

# Plot Geometry
geom.Plot()

# Plot Point Index
geom.ShowPoints()
# Plot Domain Numbers
geom.ShowDomains()

# Set Meshing Parameters
mparam = MeshingParameters()
mparam.maxh = 0.1

mesh = geom.GenerateMesh(mparam)