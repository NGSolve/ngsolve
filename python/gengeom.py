
geom = SplineGeometry()

# Define Points
pi1 = geom.AppendPoint(0,0)
pi2 = geom.AppendPoint(1,0)
pi3 = geom.AppendPoint(1,1)
pi4 = geom.AppendPoint(0,1)

# Define Segments
geom.AppendSegment([pi1,pi2])
geom.AppendSegment([pi2,pi3])
geom.AppendSegment([pi3,pi4])
geom.AppendSegment([pi4,pi1])

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