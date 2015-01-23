from ngsolve import *

pde = comp.PDE("d1_square.pde")
pde.Solve()
mesh = pde.Mesh()

lh = ngstd.LocalHeap (10000, "heap")
v = pde.spaces["v"]

lap = fem.BFI (name="laplace", dim=2, coef=fem.VariableCF ("x*y"))
src = fem.LFI (name="source", dim=2, coef=fem.ConstantCF (4.8))

i = comp.ElementId(comp.VOL,1)
el = v.GetFE(i, lh)
trafo = mesh.GetTrafo(i,lh)

mat = bla.Matrix(el.ndof, el.ndof)
lap.CalcElementMatrix (el, trafo, mat, lh)
print ("laplace matrix:\n", mat)

rhs = bla.Vector(el.ndof)
src.CalcElementVector (el, trafo, rhs, lh)
print ("source vector:\n", rhs)


