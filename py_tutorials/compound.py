from ngsolve.fem import *
from ngsolve.comp import *


mesh = Mesh("square.vol")

fes1 = FESpace("h1ho", mesh, order=3)
fes2 = FESpace("l2ho", mesh, order=2)

fes = FESpace([fes1,fes2])


f = LinearForm (fes)
f.components[0] += LFI("neumann", coef=1)
f.components[1] += LFI("source", coef=1)

f.Assemble()

print (f.vec)



a = BilinearForm (fes)
a.components[0] += BFI("mass", coef=1)
a.components[1] += BFI("mass", coef=1)

print ("integrators: ", a.integrators)
a.Assemble()
print (a.mat)







for el in fes.Elements(BND):
     print (el.GetFE().ndof)
     print (el.dofs)



