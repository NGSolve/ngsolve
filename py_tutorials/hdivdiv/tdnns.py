#
# The tangential-displacement normal-normal-stress continuous method for elasticity
#
# A.S. Pechstein and J. Schoeberl:
# Tangential-displacement and normal-normal-stress continuous mixed finite elements for elasticity.
# Mathematical Models and Methods in Applied Sciences 21(8), 1761–1782, 2011.
# 

from ngsolve import *
from netgen.geom2d import SplineGeometry

geo = SplineGeometry()
geo.AddRectangle( (0, 0), (10, 1), bcs = ("bottom", "right", "top", "left"))
mesh = Mesh( geo.GenerateMesh(maxh=0.5))

order = 3
V = HDivDiv(mesh, order=order-1, dirichlet="bottom|right|top", flags = { "plus" : True } )
Q = HCurl(mesh, order=order, dirichlet="left", flags = { "type1" : True } )
X = FESpace([V,Q])

print ("ndof-V:", V.ndof, ", ndof-Q:", Q.ndof)

sigma, u = X.TrialFunction()
tau, v = X.TestFunction()

n = specialcf.normal(2)

sigma.dims = (2,2)
tau.dims = (2,2)

def tang(u): return u-(u*n)*n

a = BilinearForm(X, symmetric=True)
a += SymbolicBFI ( InnerProduct (sigma, tau) + div(sigma)*v + div(tau)*u - 1e-10 * u*v )
a += SymbolicBFI ( -(sigma*n) * tang(v) - (tau*n)*tang(u), element_boundary=True)
a.Assemble()

f = LinearForm(X)
f += SymbolicLFI ( 1 * v[1] )
f.Assemble()

u = GridFunction(X)
u.vec.data = a.mat.Inverse(X.FreeDofs(), inverse="sparsecholesky") * f.vec

Draw (u.components[0], mesh, name="sigma")
Draw (u.components[1], mesh, name="disp")

Draw (u.components[0][0], mesh, name="s11")


