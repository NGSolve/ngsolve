from netgen.geom2d import unit_square
from ngsolve import *


ngsglobals.msg_level = 1
mesh = Mesh(unit_square.GenerateMesh(maxh=0.4))
for k in range(5):
    mesh.Refine()


order = 3
fes1 = L2(mesh, order=order)
fes2 = FacetFESpace(mesh, order=order, dirichlet=[1,2,3])

print ("element dofs: ", fes1.ndof)
print ("facet dofs: ", fes2.ndof)

fes = FESpace([fes1,fes2])

u,uhat = fes.TrialFunction()
v,vhat = fes.TestFunction()

grad_u = u.Deriv()
grad_v = v.Deriv()
n = specialcf.normal(mesh.dim)
h = specialcf.mesh_size

a = BilinearForm(fes, symmetric=True, flags = { "eliminate_internal" : True })
a += SymbolicBFI(grad(u) * grad(v))
a += SymbolicBFI(grad(u)*n*(vhat-v)+grad(v)*n*(uhat-u)+10*order*order/h*(u-uhat)*(v-vhat), element_boundary=True)

c = Preconditioner(type="direct", bf=a, flags = { "inverse" : "sparsecholesky" } )
# c = Preconditioner(type="bddc", bf=a)

a.Assemble()
ainv = CGSolver(a.mat, c.mat)

f = LinearForm(fes)
f += SymbolicLFI(1*v)
f.Assemble()

u = GridFunction(fes)


f.vec.data += a.harmonic_extension_trans * f.vec

u.vec.data = ainv * f.vec

u.vec.data += a.harmonic_extension * u.vec
u.vec.data += a.inner_solve * f.vec

Draw (u.components[0], mesh, "sol")


