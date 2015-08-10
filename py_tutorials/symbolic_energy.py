from ngsolve.fem import *
from ngsolve.la import *
from ngsolve.comp import *
from ngsolve.solve import *
from ngsolve.utils import *

import netgen
import netgen.geom2d

# mesh = Mesh (netgen.csg.unit_cube.GenerateMesh(maxh=0.2))
mp = netgen.meshing.MeshingParameters (maxh=0.2)
mesh = Mesh (netgen.geom2d.unit_square.GenerateMesh(mp))


V = H1(mesh, order=4, dirichlet=[1,2,3,4])

u = V.TrialFunction()
gradu = u.Deriv()

a = BilinearForm (V, symmetric=False)
a += SymbolicEnergy (0.05*gradu*gradu+u*u*u*u-100*u)

u = GridFunction (V)
u.vec[:] = 0

res = u.vec.CreateVector()
w = u.vec.CreateVector()

for it in range(20):
    print ("Newton iteration", it)
    print ("energy = ", a.Energy(u.vec))

    a.Apply (u.vec, res)
    a.AssembleLinearization (u.vec)
    inv = a.mat.Inverse(V.FreeDofs())
    w.data = inv * res
    print ("w*r =", InnerProduct(w,res))

    u.vec.data -= w
    Draw (u, sd=4)
    

