from ngsolve import *
from netgen.geom2d import unit_square
from solvenlmin import *
mesh = Mesh (unit_square.GenerateMesh(maxh=0.2))
V = H1(mesh, order=4, dirichlet=[1,2,3,4])
u = V.TrialFunction()

a = BilinearForm (V, symmetric=False)
a += SymbolicEnergy ( grad(u)*grad(u) + u*u*u*u-u )

gfu = GridFunction (V)
gfu.vec[:] = 0
Draw(gfu,mesh,"u")

SolveNonlinearMinProblem(a,gfu)

print ("energy = ", a.Energy(gfu.vec))
