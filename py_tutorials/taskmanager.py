from netgen.csg import unit_cube
from ngsolve import *

mesh = Mesh (unit_cube.GenerateMesh(maxh=0.4))
for l in range(3):
    mesh.Refine()

fes = H1(mesh, order=3)
u,v = fes.TnT()
a = BilinearForm(fes)
a += grad(u) * grad(v) * dx


print('sequential assembly...')
a.Assemble()

print('parallel assembly...')
with TaskManager():
    a.Assemble()
    
