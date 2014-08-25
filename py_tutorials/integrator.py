import sys
sys.path.append("/opt/netgen/lib")

from libngspy import *
from Ngsolve import *
from Ngcomp import *
from Ngfem import *
from Ngbla import *
from Ngstd import *






pde = PDE()
pde.Load ("../pde_tutorial/d1_square.pde")
pde.Solve()
mesh = pde.Mesh()

lh = LocalHeap (10000, "heap")
v = pde.spaces["v"]

lam = ConstantCF (4.8)
lap = CreateBFI (name="laplace", dim=2, coef=lam)
src = CreateLFI (name="source", dim=2, coef=lam)

i = ElementId(VOL,1)
el = v.GetFE(i, lh)
trafo = mesh.GetTrafo(i,lh)

mat = Matrix(el.ndof, el.ndof)
lap.CalcElementMatrix (el, trafo, mat, lh)
print ("laplace matrix:\n", mat)

rhs = Vector(el.ndof)
src.CalcElementVector (el, trafo, rhs, lh)
print ("source vector:\n", rhs)


