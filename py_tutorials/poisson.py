# solve the Poisson equation -Delta u = f
# with Dirichlet boundary condition u = 0

from ngsolve import *

ngsglobals.msg_level = 1

# generate a triangular mesh of mesh-size 0.2
mesh = Mesh(unit_square.GenerateMesh(maxh=0.2))

# H1-conforming finite element space
fes = H1(mesh, order=3)

# define trial- and test-functions
u,v = fes.TnT()

# the right hand side
f = 32 * (y*(1-y)+x*(1-x))

# solve the variational equation,
# Dirichlet conditions are given on named boundaries
gfu = Solve(grad(u)*grad(v)*dx == f*v*dx,
            u[BND("left|right|bottom|top")] == 0)

# plot the solution (netgen-gui only)
Draw (gfu)
Draw (-grad(gfu), mesh, "Flux")

exact = 16*x*(1-x)*y*(1-y)
print ("L2-error:", sqrt (Integrate ( (gfu-exact)*(gfu-exact), mesh)))
