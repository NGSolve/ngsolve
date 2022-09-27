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

h = 0.01
p = 5
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



print("\n*** ng ***")

n = a.mat.height

a_ngsparse = a.mat

x_ngvec = BaseVector(n * [1])

res_ngvec = BaseVector(n * [0])
t0 = time.time()
for i in range(rep):
    a_ngsparse.Mult(x_ngvec, res_ngvec)
ng_time = time.time() - t0



print("\n*** cusparse ***")

InitCuLinalg()

a_cusparse = CreateDevMatrix(a.mat)

univ = UnifiedVector(n * [1])
unires = UnifiedVector(n * [0])

univ.UpdateDevice()
unires.UpdateDevice()

t0 = time.time()
for i in range(rep):
    a_cusparse.Mult(univ, unires)
# unires.UpdateHost()
print("time uni:", time.time() - t0)


err = (res_ngvec - unires).Evaluate()
M = max([abs(err[i]) for i in range(n)])
print("max error:", M)



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
