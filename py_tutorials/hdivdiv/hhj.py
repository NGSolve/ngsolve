#
# The Hellan-Herrmann-Johnson method for a Kirchhoff plate
# 
# M. I. Comodi: The Hellan–Herrmann-Johnson Method
# Some error estimates and postprocessing. Math. Comp. 52, 17–39, 1989
# 

from ngsolve import *
from netgen.geom2d import unit_square

mesh = Mesh (unit_square.GenerateMesh(maxh=0.05))
order = 3

V = HDivDiv(mesh, order=order-1)
Q = H1(mesh, order=order, dirichlet=[1,2,3,4])
X = FESpace([V,Q])

print ("ndof-V:", V.ndof, ", ndof-Q:", Q.ndof)

sigma, u = X.TrialFunction()
tau, v = X.TestFunction()

n = specialcf.normal(2)

def tang(u): return u-(u*n)*n

a = BilinearForm(X, symmetric=True)
a += SymbolicBFI ( InnerProduct (sigma, tau) + div(sigma)*grad(v) + div(tau)*grad(u) - 1e-10*u*v )
a += SymbolicBFI ( -(sigma*n) * tang(grad(v)) - (tau*n)*tang(grad(u)), element_boundary=True)
a.Assemble()

f = LinearForm(X)
f += SymbolicLFI ( 1 * v )
# f += SymbolicLFI (tau.Trace(), BND, definedon=[1])
f.Assemble()

u = GridFunction(X)
u.vec.data = a.mat.Inverse(X.FreeDofs()) * f.vec

Draw (u.components[0], mesh, name="sigma")
Draw (u.components[1], mesh, name="disp")



