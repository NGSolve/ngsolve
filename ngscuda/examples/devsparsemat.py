# file: devsparsemat.py
# date: 28.09.2022
#
# testing basic functions for sparse matrix on device

from ngsolve import *
from ngsolve.ngscuda import *
from ngsolve.la import *
from ngsolve.webgui import Draw
from netgen.geom2d import unit_square
import time

rep = 100

h = 0.01
p = 10
mesh = Mesh(unit_square.GenerateMesh(maxh=h))
fes = H1(mesh, order=p, dirichlet="bottom|right")
fes.ndof
print("h:", h)
print("p:", p)

u, v = fes.TnT()
gfu = GridFunction(fes)

a = BilinearForm(grad(u)*grad(v)*dx, symmetric=True)
a.Assemble()

f = LinearForm(x*v*dx)
f.Assemble()


# Matrix operations on CPU
n = a.mat.height

x_cpu = a.mat.CreateVector()
x_cpu.FV()[:] = 1

res_cpu = a.mat.CreateVector()

print("\n*** Matrix-Vector Mult on the CPU ***")
a.mat.Mult(x_cpu, res_cpu)


# Matrix operations on GPU
print("\n*** Matrix-Vector Mult on the GPU ***")

a_dev = CreateDevMatrix(a.mat)

x_dev = a_dev.CreateVector()
x_dev.FV()[:] = 1

res_dev = a_dev.CreateVector()

x_dev.UpdateDevice()
res_dev.UpdateDevice()

a_dev.Mult(x_dev, res_dev)
# res_dev.UpdateHost()


# diff = (res_dev - res_cpu).Evaluate()
# M = max([abs(diff[i]) for i in range(len(diff))])
# print("Max Diff: ", M)


## timings of gpu may not be accurate due to synchronization
# print("\nTimers:")
# timers = Timers()
# for t in timers:
#     if t["name"] == "SparseMatrix::MultAdd":
#         time_cpu = t["time"]

#     if t["name"] == "CUDA Matrix-Vector Multiplication":
#         time_dev = t["time"]

# print("\nNZE / Time")
# print("NGSolve: {:.2E}".format(rep * a.mat.nze / time_cpu))
# print("Unified: {:.2E}".format(rep * a.mat.nze / time_dev))
