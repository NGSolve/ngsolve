from ngsolve import *
from math import pi

#
# solve the Poisson equation -Delta u = f
#
# with boundary conditions
#      u = 0  on Gamma1
#  du/dn = 1  on Gamma2

from netgen.geom2d import unit_square

ngsglobals.msg_level = 1
mesh = Mesh (unit_square.GenerateMesh(maxh=0.1))


# define a finite element space
# Dirichlet boundary is Gamma_1,2,3
# play around with -order=...
v = H1(mesh, order=3, dirichlet=[1,2,3])

# the solution field
u = GridFunction (space=v)

f = LinearForm (space=v)
f += Source (x*sin(pi*y))
# neumann 1 --definedon=[2]    # Neumann on Gamma_2

# the bilinear-form 
a = BilinearForm (space=v, symmetric=True) #  -eliminate_internal
a += Laplace (1)


# c = Preconditioner(bf=a, type="multigrid", { "smoother" : "block" }) 
# c = Preconditioner(bf=a, type="bddc") #  c -type=bddc -bilinearform=a
c = Preconditioner(bf=a, type="local")
# c = Preconditioner(bf=a, type="direct") 

a.Assemble()
f.Assemble()

# c.Update()

BVP (bf=a,lf=f,gf=u,pre=c,maxsteps=200).Do()

flux = -u.Deriv()

Draw (u)
Draw (flux, mesh, "Flux")

