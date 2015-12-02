from netgen.geom2d import unit_square
from ngsolve import *

ngsglobals.msg_level = 1

mesh = Mesh(unit_square.GenerateMesh(maxh=0.1))

order = 2
fes1 = HDiv(mesh, order=order)
fes2 = L2(mesh, order=order-1)

fes = FESpace([fes1,fes2])

sigma,u = fes.TrialFunction()
tau,v = fes.TestFunction()

a = BilinearForm(fes)
a += SymbolicBFI(sigma*tau + div(sigma)*v + div(tau)*u - 1e-10*u*v)
# (regularization needed for direct solver)
a.Assemble()

f = LinearForm(fes)
f += SymbolicLFI(-v)
f.Assemble()

u = GridFunction(fes)
u.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec

Draw (u.components[0], mesh, "flux")
Draw (u.components[1], mesh, "sol")
