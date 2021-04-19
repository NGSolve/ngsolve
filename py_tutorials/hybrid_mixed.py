from netgen.geom2d import unit_square
from ngsolve import *


ngsglobals.msg_level = 1

mesh = Mesh(unit_square.GenerateMesh(maxh=0.2))
# mesh.Refine()

order = 2
fes1 = HDiv(mesh, order=order, discontinuous=True)
fes2 = L2(mesh, order=order-1)
fes3 = FacetFESpace(mesh, order=order, dirichlet="top|bottom|right")

fes = fes1*fes2*fes3

sigma,u,uhat = fes.TrialFunction()
tau,v,vhat = fes.TestFunction()

n = specialcf.normal(mesh.dim)

a = BilinearForm(fes, symmetric=False, condense = True)
a += (sigma*tau + div(sigma)*v + div(tau)*u)*dx
a += (sigma*n*vhat+tau*n*uhat)*dx(element_boundary=True)

# c = Preconditioner(type="direct", bf=a, inverse = sparsecholesky)
# c = Preconditioner(type="direct", bf=a, inverse = pardiso)
# c = Preconditioner(type="local", bf=a)
c = Preconditioner(type="bddc", bf=a)


a.Assemble()

f = LinearForm(fes)
f += -1*v*dx
f.Assemble()

u = GridFunction(fes)

BVP(bf=a,lf=f,gf=u,pre=c).Do()
# u.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec

Draw (u.components[1], mesh, "sol")
Draw (u.components[0], mesh, "flux")
