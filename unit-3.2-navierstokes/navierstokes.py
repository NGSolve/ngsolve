import netgen.gui
from math import pi
from ngsolve import *
from netgen.geom2d import SplineGeometry

from netgen.geom2d import SplineGeometry
geo = SplineGeometry()
geo.AddRectangle( (0, 0), (2, 0.41), 
                 bcs = ("wall", "outlet", "wall", "inlet"))
geo.AddCircle ( (0.2, 0.2), r=0.05, 
               leftdomain=0, rightdomain=1, bc="cyl")
mesh = Mesh( geo.GenerateMesh(maxh=0.08))
mesh.Curve(3)

# viscosity
nu = 0.001


k = 3
V = H1(mesh,order=k, dirichlet="wall|cyl|inlet")
Q = H1(mesh,order=k-1)

X = FESpace([V,V,Q])

u = GridFunction(X)
velocity = CoefficientFunction (u.components[0:2])
Draw(velocity,mesh,"u")
Draw(u.components[2],mesh,"p")

# parabolic inflow at bc=1:
uin = 1.5*4*y*(0.41-y)/(0.41*0.41)
u.components[0].Set(uin, definedon=mesh.Boundaries("inlet"))

ux,uy,p = X.TrialFunction()
vx,vy,q = X.TestFunction()

div_u = grad(ux)[0]+grad(uy)[1]
div_v = grad(vx)[0]+grad(vy)[1]

stokes = nu*grad(ux)*grad(vx)+nu*grad(uy)*grad(vy)-div_u*q-div_v*p
        
a = BilinearForm(X)
a += SymbolicBFI(stokes)
a.Assemble()

f = LinearForm(X)   
f.Assemble()

inv_stokes = a.mat.Inverse(X.FreeDofs())

res = f.vec.CreateVector()
res.data = f.vec - a.mat*u.vec
u.vec.data += inv_stokes * res

dt = 0.001

# matrix for implicit Euler 
mstar = BilinearForm(X)
mstar += SymbolicBFI(ux*vx+uy*vy + dt*stokes)
mstar.Assemble()
inv = mstar.mat.Inverse(X.FreeDofs())

velocity = CoefficientFunction (u.components[0:2])
absvelocity = sqrt(velocity*velocity)

conv = LinearForm(X)
conv += SymbolicLFI(velocity*grad(u.components[0])*vx)
conv += SymbolicLFI(velocity*grad(u.components[1])*vy)

t = 0
tend = 0

# implicit Euler/explicit Euler splitting method:
tend += 1
while t < tend-0.5*dt:
    print ("\rt=", t, end="")

    conv.Assemble()
    res.data = a.mat * u.vec + conv.vec
    u.vec.data -= dt * inv * res    

    t = t + dt
    Redraw()
