import sys
sys.path.append("/opt/netgen/lib")

import libngspy
from Ngstd import *
from Ngbla import *
from Ngfem import *
from Ngcomp import *
from Ngsolve import *




pde = PDE("../pde_tutorial/d1_square.pde")
# SetDefaultPDE (pde)
mesh = pde.Mesh()

v = pde.spaces["v"]
v.Update (heapsize=1000000)

lh = LocalHeap (10000, "heap")

lam = ConstantCF (4.8)
lap = CreateBFI (name="laplace", dim=2, coef=lam)



for i in mesh.Elements(VOL):
    hr = HeapReset (lh)
    el = v.GetFE(i,lh)
    trafo = mesh.GetTrafo(i,lh)
    mat = Matrix(el.ndof, el.ndof)
    lap.CalcElementMatrix(el, trafo, mat, lh)
    print ("Element matrix of element ", i, ":\n", mat)
    hr = 0


