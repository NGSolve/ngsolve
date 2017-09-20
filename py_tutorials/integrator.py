from ngsolve import *
import ngsolve.utils as utils

pde = PDE("d1_square.pde")
pde.Solve()
mesh = pde.Mesh()

v = pde.spaces["v"]

lap = BFI (name="laplace", coef=utils.x*utils.y)
src = LFI (name="source", coef=4.8)

i = ElementId(VOL,1)
el = v.GetFE(i)
trafo = mesh.GetTrafo(i)

mat = lap.CalcElementMatrix (el, trafo)
print ("laplace matrix:\n", mat)

rhs = src.CalcElementVector (el, trafo)
print ("source vector:\n", rhs)


