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
# We use the master number as the copy argument to create a minion edge.
periodic.Append ( ["line", pnums[0], pnums[3]], leftdomain=0, rightdomain=1, copy=lright, bc="periodic")

from ngsolve import *
# from netgen.geom2d import unit_square

mesh = Mesh(periodic.GenerateMesh(maxh=0.05))


order=4
fes = L2(mesh, order=order, dgjumps = True)
u = fes.TrialFunction()
v = fes.TestFunction()

jump_u = u-u.Other()
jump_v = v-v.Other()
n = specialcf.normal(2)
mean_dudn = 0.5*n * (grad(u)+grad(u.Other()))
mean_dvdn = 0.5*n * (grad(v)+grad(v.Other()))

alpha = 4
h = specialcf.mesh_size
a = BilinearForm(fes)
a += SymbolicBFI(grad(u)*grad(v))
a += SymbolicBFI(alpha*order**2/h*jump_u*jump_v, skeleton=True)
a += SymbolicBFI(-mean_dudn*jump_v -mean_dvdn*jump_u, skeleton=True)
a += SymbolicBFI(alpha*order**2/h*u*v, BND, skeleton=True, definedon=mesh.Boundaries("outer"))
a += SymbolicBFI(-n*grad(u)*v-n*grad(v)*u, BND, skeleton=True, definedon=mesh.Boundaries("outer"))
a.Assemble()

f = LinearForm(fes)
f += SymbolicLFI(1*v)
f.Assemble()

gfu = GridFunction(fes, name="uDG")
gfu.vec.data = a.mat.Inverse() * f.vec
Draw (gfu)
