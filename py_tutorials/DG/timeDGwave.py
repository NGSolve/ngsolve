from netgen.geom2d import unit_square
from ngsolve import *

mesh = Mesh (unit_square.GenerateMesh(maxh=0.05))

order=4
fes1 = L2(mesh, order=order)
fes = fes1*fes1*fes1

p,ux,uy = fes.TrialFunction()
q,vx,vy = fes.TestFunction()

u0 = exp (-400 * ( (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) ))

n = specialcf.normal(2)

v = CoefficientFunction( (vx, vy) )
u = CoefficientFunction( (ux, uy) )

a1 = BilinearForm(fes)
a1 += grad(p) * v * dx
a1 += -0.5 * (p - p.Other()) * (v*n) * dx(element_boundary = True)

a2 = BilinearForm(fes)
a2 += -grad(q) * u * dx
a2 += 0.5 * (q - q.Other()) * (u*n) * dx(element_boundary = True)


u = GridFunction(fes)
u.components[0].Set(u0)

res = u.vec.CreateVector()
w = u.vec.CreateVector()

Draw (u.components[1], mesh, "ux")
Draw (u.components[2], mesh, "uy")
Draw (u.components[0], mesh, "p")
SetVisualization(min=-0.1, max=0.1, deformation=True)

tau = 1e-3
tend = 3

t = 0
nd = fes1.ndof

input ("<press enter>")

with TaskManager():
    while t < tend:
        a1.Apply (u.vec, w)
        fes1.SolveM (rho=CoefficientFunction(1), vec=w.Range(nd,2*nd))
        fes1.SolveM (rho=CoefficientFunction(1), vec=w.Range(2*nd,3*nd))
        u.vec.data -= tau * w

        a2.Apply (u.vec, w)
        fes1.SolveM (rho=CoefficientFunction(1), vec=w.Range(0,nd))
        u.vec.data -= tau * w

        t += tau
        Redraw()
    
