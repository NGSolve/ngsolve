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
V = HDivDiv(mesh, order=order-1, dirichlet="bottom|right|top", plus = True)
Q = HCurl(mesh, order=order, dirichlet="left", type1 = True)
X = V*Q

print ("ndof-V:", V.ndof, ", ndof-Q:", Q.ndof)

sigma, u = X.TrialFunction()
tau, v = X.TestFunction()

n = specialcf.normal(2)

def tang(u): return u-(u*n)*n

a = BilinearForm(X, symmetric=True)
a += (InnerProduct (sigma, tau) + div(sigma)*v + div(tau)*u - 1e-10 * u*v)*dx
a += (-(sigma*n) * tang(v) - (tau*n)*tang(u))*dx(element_boundary=True)
a.Assemble()

f = LinearForm(X)
f +=  1 * v[1] * dx
f.Assemble()

u = GridFunction(X)
u.vec.data = a.mat.Inverse(X.FreeDofs(), inverse="sparsecholesky") * f.vec

Draw (u.components[0], mesh, name="sigma")
Draw (u.components[1], mesh, name="disp")

Draw (u.components[0][0], mesh, name="s11")


