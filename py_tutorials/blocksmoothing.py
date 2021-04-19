from ngsolve import *
from netgen.geom2d import unit_square

mesh = Mesh(unit_square.GenerateMesh(maxh=0.2))

fes = H1(mesh, order=3)
u,v = fes.TrialFunction(), fes.TestFunction()

a = BilinearForm(fes)
a += (grad(u)*grad(v)+u*v)*dx
a.Assemble()

f = LinearForm(fes)
f += v*dx
f.Assemble()

blocks = [set() for x in range(mesh.nv)]
for el in fes.Elements():
    for v in el.vertices:
        blocks[v.nr] |= set(el.dofs)

# print (blocks)

smoother = a.mat.CreateBlockSmoother(blocks)

u = GridFunction(fes)
u.vec[:] = 0
res = f.vec.CreateVector()

for it in range(100):
    res.data = f.vec - a.mat*u.vec
    print ("|res| = ", Norm(res))
    smoother.Smooth(u.vec, f.vec, 1)
    

