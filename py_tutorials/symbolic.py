from ngsolve import *
from netgen.csg import unit_cube

ngsglobals.msg_level = 1
mesh = Mesh (unit_cube.GenerateMesh(maxh=0.2))

import timeit


V = H1(mesh, order=3)

u = V.TrialFunction()
v = V.TestFunction()

gradu = Grad(u)
gradv = Grad(v)
dvdx = gradv[0]
dvdy = gradv[1]
dvdz = gradv[2]

a = BilinearForm (V, symmetric=True)
a += (u*v+gradu*gradv)*dx
a.Assemble()


f = LinearForm (V)
f += x*v*dx
f.Assemble()

print (a.mat)
print (f.vec)
