from netgen.geom2d import *
from ngsolve import *
from solvenlmin import *

periodic = SplineGeometry()
pnts = [ (0,0), (1,0), (1,1), (0,1) ]
pnums = [periodic.AppendPoint(*p) for p in pnts]

lright = periodic.Append ( ["line", pnums[0], pnums[1]],bc="periodic")
btop = periodic.Append ( ["line", pnums[1], pnums[2]], bc="periodic")
periodic.Append ( ["line", pnums[3], pnums[2]], leftdomain=0, rightdomain=1, copy=lright, bc="periodic")
periodic.Append ( ["line", pnums[0], pnums[3]], leftdomain=0, rightdomain=1, copy=btop, bc="periodic")

mesh = Mesh (periodic.GenerateMesh(maxh=0.2))

V = Periodic(H1(mesh, order=4, dirichlet=[]))

u = V.TrialFunction()

eps = 4e-3
dt = 1e-1

gfu = GridFunction(V)
gfuold = GridFunction(V)

a = BilinearForm (V, symmetric=False)
a += SymbolicEnergy (eps/2*grad(u)*grad(u)
                     + ((1-u*u)*(1-u*u))
                     + 0.5/dt*(u-gfuold)*(u-gfuold))

from math import pi
gfu = GridFunction(V)
gfu.Set(sin(2*pi*x))
#gfu.Set(sin(1e8*x)) #<- essentially a random function
Draw(gfu,mesh,"u")
SetVisualization (deformation=False)
t = 0

for timestep in range(50):
    gfuold.vec.data = gfu.vec
    SolveNonlinearMinProblem(a,gfu)
    Redraw() 
    t += dt
    print("t = ", t)
