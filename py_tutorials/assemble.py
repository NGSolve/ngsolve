from ngsolve import *
from netgen.geom2d import unit_square

ngmesh = unit_square.GenerateMesh(maxh=0.2)
mesh = Mesh(ngmesh)

fes = H1(mesh, order=2)
u = fes.TrialFunction()
v = fes.TestFunction()
lap = SymbolicBFI(grad(u)*grad(v))

for el in fes.Elements():
    print ("el: ", el)
    mat = lap.CalcElementMatrix(el.GetFE(), el.GetTrafo())
    print ("Element matrix of element", ElementId(el), ":\n", mat)
    print ("dofs: ", el.dofs, "\n")


fes2 = L2(mesh, order=1)
for el1,el2 in zip (fes.Elements(), fes2.Elements()):
    print ("el1 dofs:", el1.dofs, "el2 dofs ", el2.dofs)

    
