from ngsolve import *
from netgen.geom2d import unit_square

ngsglobals.msg_level = 1

mesh = Mesh (unit_square.GenerateMesh(maxh=0.2))
V = H1(mesh, order=3, dirichlet=[1,2,3,4])
u,v = V.TnT()

f = LinearForm (V)
f += 32 * (y*(1-y)+x*(1-x)) * v * dx
a = BilinearForm (V, symmetric=False)
a += grad(u) * grad(v) * dx

a.Assemble()
f.Assemble()

u = GridFunction (V)
u.vec.data = a.mat.Inverse(V.FreeDofs(), inverse="umfpack") * f.vec

exact = 16*x*(1-x)*y*(1-y)
print ("L2-error:", sqrt (Integrate ( (u-exact)*(u-exact), mesh)))
