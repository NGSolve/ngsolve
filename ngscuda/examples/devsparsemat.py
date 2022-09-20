# file: devsparsemat.py
# date: 20.09.2022
#
# testing basic functions for sparse matrix on device

from ngsolve import *
from ngsolve.ngscuda import *
from ngsolve.la import *
from ngsolve.webgui import Draw
from netgen.geom2d import unit_square
import time

rep = 100

h = 0.05
p = 20
mesh = Mesh(unit_square.GenerateMesh(maxh=h))
fes = H1(mesh, order=p, dirichlet="bottom|right")
fes.ndof
print("h:", h)
print("p:", p)

u = fes.TrialFunction()  
v = fes.TestFunction()
gfu = GridFunction(fes)

u, v = fes.TnT()

a = BilinearForm(fes, symmetric=True)
a += grad(u)*grad(v)*dx
a.Assemble()

f = LinearForm(fes)
f += x*v*dx
f.Assemble()


a_ngsparse = a.mat
a_cusparse = CreateDevMatrix(a.mat)

n = a_ngsparse.height
v = UnifiedVector(n)
for i in range(n):
    v[i] = i

res_uniform = UnifiedVector(n)
for i in range(n):
    res_uniform[i] = 0


v.UpdateDevice()
res_uniform.UpdateDevice()


x_ngvec = BaseVector(n)
for i in range(n):
    x_ngvec[i] = v[i]

res_ngvec = BaseVector(n)
print("\n*** ng ***")
t0 = time.time()
for i in range(rep):
    a_ngsparse.Mult(x_ngvec, res_ngvec)
ng_time = time.time() - t0
print("time ng:", ng_time) 


print("\n*** cusparse ***")
t0 = time.time()
for i in range(rep):
    a_cusparse.Mult(v, res_uniform)
res_uniform.UpdateHost()
print("time uni:", time.time() - t0)

# M = 0
# for i in range(n):
#     tmp = res_ngvec[i] - res_uniform[i]
#     if (tmp > M):
#         M = tmp
# print("max diff: ", M)

print("\nTimers:")
timers = Timers()
for t in timers:
    if t['name'] == "CUDA MV":
       uni_time = t['time'] 
    if t['name'][:4] == "CUDA":
        print(t['name'], t['time'])

print("\nNZE / Time")
print("NGSolve: {:.2E}".format(rep * a.mat.nze / ng_time))
print("Unified: {:.2E}".format(rep * a.mat.nze / uni_time))
