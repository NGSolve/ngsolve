from netgen.geom2d import unit_square
from ngsolve import *


ngsglobals.msg_level = 1
mesh = Mesh(unit_square.GenerateMesh(maxh=0.4))
for k in range(5):
    mesh.Refine()


order = 3
fes1 = L2(mesh, order=order)
fes2 = FacetFESpace(mesh, order=order, dirichlet="bottom|right|top")

print ("element dofs: ", fes1.ndof)
print ("facet dofs: ", fes2.ndof)

fes = fes1*fes2

u,uhat = fes.TrialFunction()
v,vhat = fes.TestFunction()

n = specialcf.normal(mesh.dim)
h = specialcf.mesh_size

a = BilinearForm(fes, symmetric=True, condense = True)
a += grad(u) * grad(v) * dx
a += (grad(u)*n*(vhat-v)+grad(v)*n*(uhat-u)+10*order*order/h*(u-uhat)*(v-vhat))*dx(element_boundary=True)

c = Preconditioner(type="direct", bf=a, inverse = "sparsecholesky")
# c = Preconditioner(type="bddc", bf=a)

with TaskManager():
    a.Assemble()
    ainv = CGSolver(a.mat, c.mat)

f = LinearForm(fes)
f += 1*v*dx
f.Assemble()

u = GridFunction(fes)


f.vec.data += a.harmonic_extension_trans * f.vec

u.vec.data = ainv * f.vec

u.vec.data += a.harmonic_extension * u.vec
u.vec.data += a.inner_solve * f.vec

Draw (u.components[0], mesh, "sol")


