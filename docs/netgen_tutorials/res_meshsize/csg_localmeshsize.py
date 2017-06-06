
from netgen.csg import *

geo = CSGeometry()
# Set the mesh size on the sphere surface to 0.1
sphere = Sphere(Pnt(0,0,0),1).maxh(0.1)
# meshsize of the surface of the brick will not be finer than
# the volume mesh size
brick = OrthoBrick(Pnt(-2,-2,-2),Pnt(2,2,2))
# in the outer region we don't give a local mesh size -> global
# is used
geo.Add(brick-sphere)
# in the volume of the sphere we set the meshsize to 0.2
geo.Add(sphere,maxh=0.2)
# the global mesh size is set to 0.4
ngmesh = geo.GenerateMesh(maxh=0.4)

# for visualization we need a NGSolve mesh
from ngsolve import Mesh
Draw(Mesh(ngmesh))
