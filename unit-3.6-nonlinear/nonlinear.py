from netgen.geom2d import unit_square
from ngsolve import *

mesh = Mesh (unit_square.GenerateMesh(maxh=0.3))
V = H1(mesh, order=3, dirichlet=[1,2,3,4])
u = V.TrialFunction()
v = V.TestFunction()

a = BilinearForm(V)
a += SymbolicBFI( grad(u) * grad(v) + 5*u*u*v- 1 * v)

u = GridFunction(V)
Draw(u,mesh,"u")

res = u.vec.CreateVector()
du = u.vec.CreateVector()

for it in range(25):
    print ("Iteration",it)
    a.Apply(u.vec, res)
    a.AssembleLinearization(u.vec)
    
    du.data = a.mat.Inverse(V.FreeDofs()) * res
    u.vec.data -= du
    
    #stopping criteria
    stopcritval = sqrt(abs(InnerProduct(du,res)))
    print ("<A u",it,", A u",it,">_{-1}^0.5 = ", stopcritval)
    if stopcritval < 1e-13:
        break
