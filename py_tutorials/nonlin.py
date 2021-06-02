from netgen.geom2d import unit_square
from ngsolve import *

m = Mesh (unit_square.GenerateMesh(maxh=0.3))

V = H1(m, order=3, dirichlet="left|right|top|bottom")
u = V.TrialFunction()
v = V.TestFunction()

a = BilinearForm(V)
a += ( grad(u) * grad(v) + 5*u*u*v- 1 * v)*dx

u = GridFunction(V)
r = u.vec.CreateVector()
w = u.vec.CreateVector()

for it in range(5):
    print ("Iteration",it)
    a.Apply(u.vec, r)
    a.AssembleLinearization(u.vec)
    
    w.data = a.mat.Inverse(V.FreeDofs()) * r.data
    print ("|w| =", w.Norm())
    u.vec.data -= w

    Draw(u)
    input("<press a key>")
    
    
