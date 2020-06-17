from netgen.geom2d import *

periodic = SplineGeometry()
pnts = [ (0,0), (1,0), (1,1), (0,1) ]
pnums = [periodic.AppendPoint(*p) for p in pnts]

periodic.Append ( ["line", pnums[0], pnums[1]],bc="outer")
# This should be our master edge so we need to save its number.
lright = periodic.Append ( ["line", pnums[1], pnums[2]], bc="periodic")
periodic.Append ( ["line", pnums[2], pnums[3]], bc="outer")
# Minion boundaries must be defined in the same direction as master ones,
# this is why the the point numbers of this spline are defined in the reverse direction,
# leftdomain and rightdomain must therefore be switched as well!
# We use the master number as the copy argument to create a slave edge.
periodic.Append ( ["line", pnums[0], pnums[3]], leftdomain=0, rightdomain=1, copy=lright, bc="periodic")

from ngsolve import *
mesh = Mesh(periodic.GenerateMesh(maxh=0.2))

fes = Periodic(H1(mesh,order=3,dirichlet="outer"))

u,v = fes.TrialFunction(), fes.TestFunction()

a = BilinearForm(fes,symmetric=True)
a += SymbolicBFI(grad(u) * grad(v))

f = LinearForm(fes)
f += SymbolicLFI(x*v)

u = GridFunction(fes,"u")

with TaskManager():
    a.Assemble()
    f.Assemble()
    u.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec

Draw (u)

