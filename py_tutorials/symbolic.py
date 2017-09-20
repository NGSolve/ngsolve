from ngsolve import *
from netgen.csg import unit_cube

ngsglobals.msg_level = 1
mesh = Mesh (unit_cube.GenerateMesh(maxh=0.2))

import timeit


V = H1(mesh, order=3)

u = V.TrialFunction()
v = V.TestFunction()

gradu = u.Deriv()
gradv = v.Deriv()
dvdx = gradv[0]
dvdy = gradv[1]
dvdz = gradv[2]

a = BilinearForm (V, symmetric=True)
a += SymbolicBFI (u*v+gradu*gradv)
a.Assemble(heapsize=10000000)


f = LinearForm (V)
f += SymbolicLFI (x*v)
f.Assemble()

print (a.mat)
print (f.vec)
