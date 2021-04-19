from ngsolve import *
from netgen.geom2d import unit_square

mesh = Mesh (unit_square.GenerateMesh(maxh=0.2))


V = H1(mesh, order=4, dirichlet="left|right|top|bottom")

u = V.TrialFunction()

a = BilinearForm (V, symmetric=False)
a += Variation( (0.05*grad(u)*grad(u) + u*u*u*u - 100*u)*dx )

u = GridFunction (V)
u.vec[:] = 0

res = u.vec.CreateVector()
w = u.vec.CreateVector()

Draw(u,sd=4)

for it in range(20):
    print ("Newton iteration", it)
    print ("energy = ", a.Energy(u.vec))

    a.Apply (u.vec, res)
    a.AssembleLinearization (u.vec)
    inv = a.mat.Inverse(V.FreeDofs())
    w.data = inv * res
    print ("w*r =", InnerProduct(w,res))

    u.vec.data -= w
    Redraw()
    

