from netgen.geom2d import unit_square
from ngsolve import *
from math import pi
import scipy.linalg
from scipy import random
from random import random

mesh = Mesh(unit_square.GenerateMesh(maxh=0.2))

fes = H1(mesh, order=8, dirichlet=[1,2,3,4])
u = fes.TrialFunction()
v = fes.TestFunction()

a = BilinearForm(fes)
a += SymbolicBFI(grad(u)*grad(v))
pre = Preconditioner(a, "multigrid")

m = BilinearForm(fes)
m += SymbolicBFI(u*v)

a.Assemble()
m.Assemble()

num = 8
u = GridFunction(fes, multidim=num)

r = u.vec.CreateVector()
Av = u.vec.CreateVector()
Mv = u.vec.CreateVector()

vecs = []
for i in range(2*num):
    vecs.append (u.vec.CreateVector())

freedofs = fes.FreeDofs()    
for v in u.vecs:
    for i in range(len(u.vec)):
        v[i] = random() if freedofs[i] else 0

lams = num * [1]

asmall = Matrix(2*num, 2*num)
msmall = Matrix(2*num, 2*num)

for i in range(100):
    
    for j in range(num):
        vecs[j].data = u.vecs[j]
        r.data = a.mat * vecs[j] - lams[j] * m.mat * vecs[j]
        # r.data = 1/Norm(r) * r
        r *= 1/Norm(r)
        vecs[num+j].data = pre.mat * r

    for j in range(2*num):
        Av.data = a.mat * vecs[j]
        Mv.data = m.mat * vecs[j]
        for k in range(2*num):
            asmall[j,k] = InnerProduct(Av, vecs[k])
            msmall[j,k] = InnerProduct(Mv, vecs[k])

    ev,evec = scipy.linalg.eigh(a=asmall, b=msmall)
    lams[0:num] = ev[0:num]
    print (i, ":", [lam/pi**2 for lam in lams])
    
    for j in range(num):
        r[:] = 0.0
        for k in range(2*num):
            r.data += float(evec[k,j]) * vecs[k]
        u.vecs[j].data = r
    
Draw (u)
