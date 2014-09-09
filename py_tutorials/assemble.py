import sys
import os
from os import environ
sys.path.append(environ['NETGENDIR']+"/../lib")

from libngspy import *

# being sloppy ....
# from libngspy.ngstd import *
# from libngspy.ngbla import *
# from libngspy.ngfem import *
# from libngspy.ngcomp import *
# from libngspy.ngsolve import *



pde = ngsolve.PDE("../pde_tutorial/d1_square.pde")
# SetDefaultPDE (pde)
mesh = pde.Mesh()

v = pde.spaces["v"]
v.Update (heapsize=1000000)

lh = ngstd.LocalHeap (10000, "heap")

lam = ngfem.ConstantCF (4.8)
lap = ngfem.CreateBFI (name="laplace", dim=2, coef=lam)



for i in mesh.Elements():
    hr = ngstd.HeapReset(lh)
    el = v.GetFE(i,lh)
    trafo = mesh.GetTrafo(i,lh)
    mat = ngbla.Matrix(el.ndof, el.ndof)
    lap.CalcElementMatrix(el, trafo, mat, lh)
    print ("Element matrix of element ", i, ":\n", mat)
    del hr


