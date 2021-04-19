from netgen.geom2d import unit_square
from ngsolve import *

mesh = Mesh (unit_square.GenerateMesh(maxh=0.05))

order=2
fes = L2(mesh, order=order)

u = fes.TrialFunction()
v = fes.TestFunction()

u0 = exp (-40 * ( (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) ))

n = specialcf.normal(2)
h = specialcf.mesh_size

a = BilinearForm(fes)
a += grad(u) * grad(v) * dx
cf1 = -0.5 * InnerProduct(grad(u), n)*(v-v.Other(bnd=0))
cf2 = -0.5 * InnerProduct(grad(v), n)*(u-u.Other(bnd=u0))
cf3 = 2*( (order+1)**2)/h * (u-u.Other(bnd=u0)) * v
a += (cf1+cf2+cf3)*dx(element_boundary=True)

u = GridFunction(fes)
u.Set(u0)

w = u.vec.CreateVector()

Draw (u, mesh, "u")

tau = 2e-6
tend = 0.5

t = 0
with TaskManager():
    while t < tend:
        a.Apply (u.vec, w)
        fes.SolveM (rho=CoefficientFunction(1), vec=w)
        u.vec.data -= tau * w
        t += tau
        Redraw()
    
