from netgen.geom2d import unit_square
from ngsolve import *


ngsglobals.msg_level = 1
mesh = Mesh(unit_square.GenerateMesh(maxh=0.4))
mesh.Refine()
mesh.Refine()
# mesh.Refine()


order = 3
fes1 = L2(mesh, order=order)
fes2 = FESpace("facet", mesh, order=order, dirichlet=[1,2,3])

fes = FESpace([fes1,fes2])

u,uhat = fes.TrialFunction()
v,vhat = fes.TestFunction()

grad_u = u.Deriv()
grad_v = v.Deriv()
n = specialcf.normal(mesh.dim)
h = specialcf.mesh_size

a = BilinearForm(fes, symmetric=True, flags = { "eliminate_internal" : True, "printelmat" : False, "print" : False }) 
a += SymbolicBFI(grad_u * grad_v)
a += SymbolicBFI(grad_u*n*(vhat-v)+grad_v*n*(uhat-u)+10*order*order/h*(u-uhat)*(v-vhat), element_boundary=True)

c = Preconditioner(type="direct", bf=a, flags = { "inverse" : "sparsecholesky" } )
# c = Preconditioner(type="direct", bf=a, flags = { "inverse" : "pardiso" } )
# c = Preconditioner(type="local", bf=a)
# c = Preconditioner(type="bddc", bf=a)

a.Assemble()
# print (a.mat)

f = LinearForm(fes)
f += SymbolicLFI(1*v)
f.Assemble()

u = GridFunction(fes)

BVP(bf=a,lf=f,gf=u,pre=c,maxsteps=1000).Do()
# u.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec

Draw (u.components[0], mesh, "sol")


