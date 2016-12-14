from ngsolve import *
from netgen.geom2d import unit_square

mesh = Mesh(unit_square.GenerateMesh(maxh=0.05))

order = 3

V = H1(mesh, order=order)
u = V.TrialFunction()
v = V.TestFunction()

n = specialcf.normal(2)
h = specialcf.mesh_size
penalty = 3*order**2

udir = CoefficientFunction(x*y)

a = BilinearForm(V)
a += SymbolicBFI (grad(u)*grad(v))
a += SymbolicBFI (-n*grad(u)*v - n*grad(v)*u + penalty/h*u*v, BND, skeleton=True)
a.Assemble()

f = LinearForm(V)
f += SymbolicLFI ( 1 * v)
f += SymbolicLFI ( -n*grad(v)*udir + penalty/h*udir*v, BND, skeleton=True)
f.Assemble()

u = GridFunction(V)
u.vec.data = a.mat.Inverse() * f.vec

Draw (u)


