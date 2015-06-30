##import sys
##import os
##from os import environ
##sys.path.append(environ['NETGENDIR']+"/../lib")

##from libngspy import *

# being sloppy ....
# from libngspy.ngstd import *
# from libngspy.ngbla import *
# from libngspy.ngfem import *
# from libngspy.ngcomp import *
# from libngspy.ngsolve import *


from ngsolve import *

pde = comp.PDE("d1_square.pde")
# SetDefaultPDE (pde)
mesh = pde.Mesh()

v = pde.spaces["v"]
v.Update (heapsize=10000)

lap = fem.BFI ("laplace", coef=4.8)

for el in v.Elements():
    print ("el: ", el)
    mat = lap.CalcElementMatrix(el.GetFE(), el.GetTrafo())
    print ("Element matrix of element ", el, ":\n", mat)
    print ("dofs: ", el.dofs, "\n")



v2 = comp.FESpace ("h1ho", mesh, order=1, dirichlet=[1])
v2.Update()

for el1,el2 in zip (v.Elements(), v2.Elements()):
    print ("el1 dofs:", el1.dofs, "el2 dofs ", el2.dofs)

