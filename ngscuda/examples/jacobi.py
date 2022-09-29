# file: jacobi.py
# date: 28.09.2022
# 
# CG-Solver using Jacobi on GPU

from ngsolve import *
import ngsolve
from ngsolve.la import *
from ngsolve.ngscuda import *
from ngsolve.webgui import Draw
from netgen.geom2d import unit_square
import time

h = 0.01
p = 5
rep = 3
mesh = Mesh(unit_square.GenerateMesh(maxh=h))
fes = H1(mesh, order=p, dirichlet="bottom|right")
fes.ndof

print("h:", h)
print("p:", p)

u, v = fes.TnT()

a = BilinearForm(grad(u)*grad(v)*dx)
a.Assemble()

f = LinearForm(x*v*dx)
f.Assemble()

gfu = GridFunction(fes)


# CG-Solver on CPU
c = a.mat.CreateSmoother(fes.FreeDofs())

with TaskManager():
    inv_host = CGSolver(a.mat, c, maxsteps=1000)
    gfu.vec.data = inv_host * f.vec


# CG-Solver on GPU
fdev = UnifiedVector(f.vec)
adev = CreateDevMatrix(a.mat)

devjac = CreateDevSmoother(a.mat, fes.FreeDofs())

fdev.UpdateDevice()

inv = CGSolver(adev, devjac, maxsteps=1000)
res = inv * fdev
res.Evaluate()


# diff = (res - gfu.vec).Evaluate()
# M = max([abs(diff[i]) for i in range(len(diff))])
# print("max diff:", M)

