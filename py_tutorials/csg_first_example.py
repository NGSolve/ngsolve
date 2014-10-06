from netgen.csg import *
import netgen.meshing as meshing


sp1 = Sphere (Point3d(1,0,0), 1)
sp2 = Sphere (Point3d(2,0,0), 1)

both = Or(sp1, sp2)
all = And (both, OrthoBrick ( Point3d(0,0,0), Point3d (3,0.5,1)))


geom = CSGeometry()
geom.Add (all)


vis = VS(geom)
vis.Draw()


# param = meshing.MeshingParameters(maxh=0.2)
# mesh = GenerateMesh (geom, param)
# mesh.Save ("test.vol")
