from ngsolve import *
from netgen.occ import *

# viscosity
nu = 0.001

# timestepping parameters
tau = 0.001
tend = 10

# channel with a cylinder, boundaries get names
channel = Rectangle(2, 0.41).Face()
channel.edges.name = "wall"
channel.edges.Min(X).name = "inlet"
channel.edges.Max(X).name = "outlet"
cyl = Circle((0.2,0.2), 0.05).Face()
cyl.edges.name = "cyl"
cyl.edges.maxh = 0.02
mesh = (channel-cyl).GenerateMesh(dim=2, maxh=0.07)

mesh.Curve(3)

V = VectorH1(mesh,order=3, dirichlet="wall|cyl|inlet")
Q = H1(mesh,order=2)

fes = V*Q

u,p = fes.TrialFunction()
v,q = fes.TestFunction()

stokes = (nu*InnerProduct(grad(u), grad(v))+div(u)*q+div(v)*p - 1e-10*p*q)*dx
a = BilinearForm(stokes).Assemble()

# parabolic inflow at inlet:
uin = CF( (1.5*4*y*(0.41-y)/(0.41*0.41), 0) )

# solve Stokes problem for initial conditions:
gfu = Solve(stokes == 0,
            u[BND("wall|cyl")] == (0,0),
            u[BND("inlet")] == uin)

# matrix for implicit Euler 
mstar = BilinearForm(u*v*dx + tau*stokes).Assemble()
inv = mstar.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")

# the non-linear term 
conv = BilinearForm((grad(u) * u) * v * dx, nonassemble = True)

# for visualization
Draw (Norm(gfu.components[0]), mesh, "velocity", sd=3)

# implicit Euler/explicit Euler splitting method:
t = 0
with TaskManager():
    while t < tend:
        print ("t=", t, end="\r")

        res = conv.Apply(gfu.vec) + a.mat*gfu.vec
        gfu.vec.data -= tau * inv * res    

        t = t + tau
        Redraw()
