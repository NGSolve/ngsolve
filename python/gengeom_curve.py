
geom = SplineGeometry()

# Define Points
pi1 = geom.AppendPoint(0,0)
pi2 = geom.AppendPoint(1,0)
pi3 = geom.AppendPoint(1,0.5)
pi4 = geom.AppendPoint(1,1)
pi5 = geom.AppendPoint(0.5,1)
pi6 = geom.AppendPoint(0,1)

# Define Segments
geom.Append(Line(pi1,pi2))
geom.Append(Line(pi2,pi3))
geom.Append(Spline3(pi3,pi4,pi5))
geom.Append(Line(pi5,pi6))
geom.Append(Line(pi6,pi1))

# Plot Geometry
geom.Plot()

# Plot Point Index
geom.ShowPoints()
# Plot Domain Numbers
geom.ShowDomains()

# Hide point indices and domain numbers
geom.ShowPoints(False)
geom.ShowDomains(False)

# Set Meshing Parameters
mparam = MeshingParameters()
mparam.maxh = 0.1

mesh = geom.GenerateMesh(mparam)