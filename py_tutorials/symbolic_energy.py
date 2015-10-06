from ngsolve import *
from netgen.geom2d import unit_square

mesh = Mesh (unit_square.GenerateMesh(maxh=0.2))


V = H1(mesh, order=4, dirichlet=[1,2,3,4])

u = V.TrialFunction()
gradu = u.Deriv()

a = BilinearForm (V, symmetric=False)
a += SymbolicEnergy (0.05*gradu*gradu+u*u*u*u-100*u)

u = GridFunction (V)
u.vec[:] = 0

res = u.vec.CreateVector()
w = u.vec.CreateVector()

for it in range(20):
    print ("Newton iteration", it)
    print ("energy = ", a.Energy(u.vec))

    a.Apply (u.vec, res)
    a.AssembleLinearization (u.vec)
    inv = a.mat.Inverse(V.FreeDofs())
    w.data = inv * res
    print ("w*r =", InnerProduct(w,res))

    u.vec.data -= w
    Draw (u, sd=4)
    

