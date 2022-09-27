from ngsolve import *
import ngsolve
from ngsolve.la import *
from ngsolve.ngscuda import *
from ngsolve.webgui import Draw
from netgen.geom2d import unit_square
import time

h = 0.01
p = 2
rep = 10
mesh = Mesh(unit_square.GenerateMesh(maxh=h))
fes = H1(mesh, order=p, dirichlet="bottom|right")
fes.ndof
print("h:", h)
print("p:", p)

u, v = fes.TnT()

a = BilinearForm(grad(u)*grad(v)*dx, symmetric=True)
# c = Preconditioner(a, "local")

a.Assemble()
c = a.mat.CreateSmoother(fes.FreeDofs())

f = LinearForm(x*v*dx)
f.Assemble()

gfu = GridFunction(fes)

t = time.time()
for r in range(rep):
    print("rep: ", r)
    inv_host = CGSolver(a.mat, c, maxsteps=1000)
    gfu.vec.data = inv_host * f.vec
print("CG time:", time.time() - t)



InitCuLinalg()

fdev = UnifiedVector(f.vec)
adev = CreateDevMatrix(a.mat)

devjac = CreateDevSmoother(a.mat, fes.FreeDofs())

fdev.UpdateDevice()

t = time.time()
for r in range(rep):
    print("rep: ", r)
    inv = CGSolver(adev, devjac, maxsteps=1000)
    res = inv * fdev
    res.Evaluate()
print("dev cg:", time.time() -t)

n = a.mat.height
err = (res - gfu.vec).Evaluate()
M = max([err[i] for i in range(n)])
print("max err:", M)
