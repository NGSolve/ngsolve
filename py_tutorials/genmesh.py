from ngsolve.ngstd import *
from ngsolve.comp import *
from ngsolve.solve import *
from ngsolve.utils import *

ngsglobals.msg_level = 1

# a pre-defined csg geometry:
from netgen.csg import unit_cube

# generate a netgen-mesh
ngmesh = unit_cube.GenerateMesh (maxh=0.2)


# convert to ngsolve-mesh
mesh = Mesh(ngmesh)

v = H1(mesh, order=3)
u = GridFunction(v)
u.Set (x*y*z)

Draw (u)
print ("Integral(u) = ", Integrate(u, mesh, order=3))







