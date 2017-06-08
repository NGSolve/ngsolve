from netgen.geom2d import unit_square
from ngsolve import *
from math import pi
import scipy.linalg
from scipy import random 
from random import random

mesh = Mesh(unit_square.GenerateMesh(maxh=0.1))

# setup stiffness and mass matrix:
fes = H1(mesh, order=5, dirichlet=[1,2,3,4])
u = fes.TrialFunction()
v = fes.TestFunction()

a = BilinearForm(fes)
a += SymbolicBFI(grad(u)*grad(v))
pre = Preconditioner(a, "multigrid")

m = BilinearForm(fes)
m += SymbolicBFI(u*v)

a.Assemble()
m.Assemble()

u = GridFunction(fes)

# some help vectors
r = u.vec.CreateVector()
w = u.vec.CreateVector()
Mu = u.vec.CreateVector()
Au = u.vec.CreateVector()
Mw = u.vec.CreateVector()
Aw = u.vec.CreateVector()

# random initial conditions, zero on boundary
freedofs = fes.FreeDofs()    
for i in range(len(u.vec)):
    u.vec[i] = random() if freedofs[i] else 0    


for i in range(100):
    Au.data = a.mat * u.vec
    Mu.data = m.mat * u.vec
    auu = InnerProduct(Au, u.vec)
    muu = InnerProduct(Mu, u.vec)
    # Rayleigh quotiont
    lam = auu/muu
    print (lam / (pi**2))
    # residual
    r.data = Au - lam * Mu
    w.data = pre.mat * r.data
    w.data = 1/Norm(w) * w
    Aw.data = a.mat * w
    Mw.data = m.mat * w

    # setup and solve 2x2 small eigenvalue problem
    asmall = Matrix(2,2)
    asmall[0,0] = auu
    asmall[0,1] = asmall[1,0] = InnerProduct(Au, w)
    asmall[1,1] = InnerProduct(Aw, w)
    msmall = Matrix(2,2)
    msmall[0,0] = muu
    msmall[0,1] = msmall[1,0] = InnerProduct(Mu, w)
    msmall[1,1] = InnerProduct(Mw, w)
    # print ("asmall =", asmall, ", msmall = ", msmall)
    eval,evec = scipy.linalg.eigh(asmall, b=msmall)
    # print (eval / pi**2)
    print (evec)
    u.vec.data = float(evec[0,0]) * u.vec + float(evec[1,0]) * w
    
Draw (u)
