import netgen.gui
import tkinter
from math import pi
from ngsolve import *
from netgen.geom2d import SplineGeometry

geo = SplineGeometry()
geo.AddRectangle( (-1, -1), (1, 1), bcs = ("bottom", "right", "top", "left"))
mesh = Mesh( geo.GenerateMesh(maxh=0.25))
fes = H1(mesh, order=3, dirichlet="bottom|right|left|top")
u = fes.TrialFunction()
v = fes.TestFunction()

time = 0.0
dt = 0.001

b = CoefficientFunction((2*y*(1-x*x),-2*x*(1-y*y)))

a = BilinearForm(fes, symmetric=False)
a += SymbolicBFI (0.01*grad(u)*grad(v) + b*grad(u)*v)
a.Assemble()

m = BilinearForm(fes, symmetric=False)
m += SymbolicBFI (u*v)
m.Assemble()

mstar = m.mat.CreateMatrix()
mstar.AsVector().data = m.mat.AsVector() + dt * a.mat.AsVector()
invmstar = mstar.Inverse(freedofs=fes.FreeDofs())

f = LinearForm(fes)
gaussp = exp(-6*((x+0.5)*(x+0.5)+y*y))-exp(-6*((x-0.5)*(x-0.5)+y*y))
f += SymbolicLFI(gaussp*v)
f.Assemble()

gfu = GridFunction(fes)
gfu.Set((1-y*y)*x)
Draw(gfu,mesh,"u")

tend = 10 
res = gfu.vec.CreateVector()
while time < tend - 0.5 * dt:
    res.data = dt * f.vec - dt * a.mat * gfu.vec
    gfu.vec.data += invmstar * res
    time += dt
    print("\r",time,end="")
    Redraw(blocking=False)
