from ngsolve import *

pde = PDE("d1_square.pde")
mesh = pde.Mesh()

v = pde.spaces["v"]
v.Update (heapsize=10000)

lap = BFI ("laplace", coef=4.8)

for el in v.Elements():
    print ("el: ", el)
    mat = lap.CalcElementMatrix(el.GetFE(), el.GetTrafo())
    print ("Element matrix of element", ElementId(el), ":\n", mat)
    print ("dofs: ", el.dofs, "\n")



v2 = H1(mesh, order=1, dirichlet=[1])
v2.Update()

for el1,el2 in zip (v.Elements(), v2.Elements()):
    print ("el1 dofs:", el1.dofs, "el2 dofs ", el2.dofs)

