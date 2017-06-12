from math import pi
from ngsolve import *
from netgen.geom2d import SplineGeometry
from netgen.geom2d import SplineGeometry
geo = SplineGeometry()
geo.AddRectangle( (0, 0), (1, 1), 
                 bcs = ("bottom", "right", "top", "left"))
mesh = Mesh( geo.GenerateMesh(maxh=0.2))

b = CoefficientFunction((1,2))
n = specialcf.normal(mesh.dim)
ubnd = IfPos(x,1,0)

V = L2(mesh,order=2)
u,v = V.TrialFunction(), V.TestFunction()

c = BilinearForm(V)
c += SymbolicBFI( b * grad(u) * v)
c += SymbolicBFI( IfPos( (b*n), 0, (b*n) * (u.Other(ubnd)-u)) * v, element_boundary=True)

gfu_expl = GridFunction(V)
Draw(gfu_expl,mesh,"u_explicit")

res = gfu_expl.vec.CreateVector()
c.Apply(gfu_expl.vec,res)

t = 0
dt = 0.001
tend = 1

while t < tend-0.5*dt:
    c.Apply(gfu_expl.vec,res)
    V.SolveM(CoefficientFunction(1.0),res)
    gfu_expl.vec.data -= dt * res
    t += dt
    Redraw(blocking=False)
