import sys
import os
from os import environ
sys.path.append(environ['NETGENDIR']+"/../lib")

from libngspy import *

pde = ngsolve.PDE("../pde_tutorial/d1_square.pde")
pde.Solve()
mesh = pde.Mesh()

lh = ngstd.LocalHeap (10000, "heap")
v = pde.spaces["v"]

lap = ngfem.BFI (name="laplace", dim=2, coef=ngfem.VariableCF ("x*y"))
src = ngfem.LFI (name="source", dim=2, coef=ngfem.ConstantCF (4.8))

i = ngcomp.ElementId(ngcomp.VOL,1)
el = v.GetFE(i, lh)
trafo = mesh.GetTrafo(i,lh)

mat = ngbla.Matrix(el.ndof, el.ndof)
lap.CalcElementMatrix (el, trafo, mat, lh)
print ("laplace matrix:\n", mat)

rhs = ngbla.Vector(el.ndof)
src.CalcElementVector (el, trafo, rhs, lh)
print ("source vector:\n", rhs)


