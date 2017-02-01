from netgen.geom2d import unit_square
from ngsolve import *


ngsglobals.msg_level = 1

mesh = Mesh(unit_square.GenerateMesh(maxh=0.2))
# mesh.Refine()

order = 2
fes1 = HDiv(mesh, order=order, flags = { "discontinuous" : True } )
fes2 = L2(mesh, order=order-1)
fes3 = FacetFESpace(mesh, order=order, dirichlet=[1,2,3])

fes = FESpace([fes1,fes2,fes3])

sigma,u,uhat = fes.TrialFunction()
tau,v,vhat = fes.TestFunction()

n = specialcf.normal(mesh.dim)

a = BilinearForm(fes, symmetric=False, flags = { "eliminate_internal" : True }) 
a += SymbolicBFI(sigma*tau + div(sigma)*v + div(tau)*u)
a += SymbolicBFI(sigma*n*vhat+tau*n*uhat, element_boundary=True)

# c = Preconditioner(type="direct", bf=a, flags = { "inverse" : "sparsecholesky" } )
# c = Preconditioner(type="direct", bf=a, flags = { "inverse" : "pardiso" } )
# c = Preconditioner(type="local", bf=a)
c = Preconditioner(type="bddc", bf=a)


a.Assemble()

f = LinearForm(fes)
f += SymbolicLFI(-1*v)
f.Assemble()

u = GridFunction(fes)

BVP(bf=a,lf=f,gf=u,pre=c).Do()
# u.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec

Draw (u.components[1], mesh, "sol")
Draw (u.components[0], mesh, "flux")
